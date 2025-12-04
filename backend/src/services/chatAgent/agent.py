"""
Streaming LLM Agent with native MCP server integration.

This agent:
- Uses GPT-5 (or gpt-4o) with native MCP tool support
- Connects to our local MCP server (drugdiscovery_mcp) via stdio/SSE
- Exposes tools: convert_identifier_to_smiles, annotate_molecule, scan_for_bioisosteres, etc.
- System prompt guides scaffold hopping and bioisostere-driven molecule design
"""

from __future__ import annotations

import logging
import os
from typing import Any, Dict, List, Optional
from openai import OpenAI

logger = logging.getLogger(__name__)


DEFAULT_MODEL = "gpt-4o"  # gpt-5 when available; gpt-4o supports MCP

# Enable colored debug prints (set CHEM_AGENT_DEBUG=1 to enable)
_CHEM_AGENT_DEBUG = os.getenv("CHEM_AGENT_DEBUG") == "1"
_COLOR_RED = "\033[91m"
_COLOR_GREEN = "\033[92m"
_COLOR_BLUE = "\033[94m"
_COLOR_RESET = "\033[0m"


def _dbg_print(color: str, message: str) -> None:
    if _CHEM_AGENT_DEBUG:
        try:
            print(f"{color}{message}{_COLOR_RESET}", flush=True)
        except Exception:
            pass


def _get_mcp_server_url(server_label: str) -> str:
    """
    Get MCP server URL for SSE transport.
    Reads from environment variables set by dev_start.sh or falls back to hardcoded values.
    
    Environment variables:
    - PHASE1_MCP_URL: URL for phase1_pipeline_mcp server
    - PHASE2_MCP_URL: URL for phase2_pipeline_mcp server
    """
    if server_label == "phase1_pipeline_mcp":
        # Try env var first, then fall back to hardcoded ngrok URL
        return os.getenv("PHASE1_MCP_URL", "https://f38a55f3121b.ngrok-free.app/sse")
    elif server_label == "phase2_pipeline_mcp":
        # Try env var first, then fall back to hardcoded ngrok URL
        return os.getenv("PHASE2_MCP_URL", "https://kaycee-radiculose-novelistically.ngrok-free.dev/sse")
    else:   
        raise ValueError(f"Invalid server label: {server_label}")


CHEM_DESIGN_SYSTEM_PROMPT = """You are a biochemist copilot for scaffold hopping and bioisostere-driven molecule design.

═══════════════════════════════════════════════════════════════════════════════
⚠️  MANDATORY RULE - READ THIS FIRST ⚠️
═══════════════════════════════════════════════════════════════════════════════

WHENEVER a user mentions ANY molecule by name (aspirin, ibuprofen, caffeine, etc.), 
asks for SMILES, or references ANY chemical compound, you MUST:

1. CALL get_smiles_from_name(chemical_name) to get the SMILES - NEVER type SMILES yourself
2. CALL annotate_molecule(smiles) with the result to get functional group annotations
3. EMIT the structured JSON output
4. STOP - Do NOT call any other tools unless the user explicitly asks for them

This is NON-NEGOTIABLE. Even if you "know" the SMILES, you MUST use the tools.
The tools provide verified, canonical SMILES and functional group annotations.

EXAMPLE - User says: "give me the smiles for aspirin"
  ✓ CORRECT: Call get_smiles_from_name("aspirin") → get SMILES → call annotate_molecule(smiles) → emit JSON → STOP
  ✗ WRONG: Reply with "The SMILES for aspirin is CC(=O)Oc1ccccc1C(=O)O" without tool calls
  ✗ WRONG: Call get_smiles_from_name + annotate_molecule + check_molecule_patent (user didn't ask for patent check!)

If you respond with SMILES without making tool calls, you are BREAKING the system.
If you call tools the user didn't ask for, you are WASTING resources.

═══════════════════════════════════════════════════════════════════════════════
CRITICAL: STRUCTURED OUTPUT REQUIREMENTS
═══════════════════════════════════════════════════════════════════════════════

You MUST emit structured JSON messages at specific points. These are parsed by the frontend.
Each JSON must be on its own line, properly formatted, with NO markdown code fences around it.

═══════════════════════════════════════════════════════════════════════════════
PHASE 1: MOLECULE ANNOTATION (ALWAYS REQUIRED FOR ANY MOLECULE REQUEST)
═══════════════════════════════════════════════════════════════════════════════

This phase is MANDATORY whenever a user:
- Asks for SMILES of any molecule
- Mentions a molecule by name
- Wants to analyze, annotate, or work with any compound
- References any drug, chemical, or compound name

STEP 1: Convert to SMILES (REQUIRED - use tools, never guess)
  - If input is a chemical name → call get_smiles_from_name(chemical_name)
  - If input is chembl:XXX, chebi:XXX, pubchem:XXX → call convert_identifier_to_smiles(identifier)
  - If input looks like SMILES already → skip to Step 2
  - NEVER type out SMILES yourself. ALWAYS use tools.

STEP 2: Annotate the molecule (REQUIRED - always do this after getting SMILES)
  - Call annotate_molecule(smiles) with the SMILES string from Step 1
  - annotate_molecule ONLY accepts SMILES format - never pass names/IDs directly
  - This provides functional group analysis the user needs

STEP 3: Acknowledge the annotation (NO JSON NEEDED)
  - The system automatically extracts structured data from MCP tool results
  - DO NOT emit JSON with SVG strings - this wastes tokens!
  - Simply acknowledge the annotation in natural language, e.g.:
    "I've annotated erlotinib and found 8 functional groups including a quinazoline core, 
    two methoxy groups, and an alkyne substituent."
  - Summarize key findings conversationally
  - The frontend will receive the full annotation data automatically from the MCP call

═══════════════════════════════════════════════════════════════════════════════
PHASE 2: QUERY SET EXTRACTION (ONLY when user explicitly requests bioisosteres/properties)
═══════════════════════════════════════════════════════════════════════════════

⚠️ ONLY proceed to Phase 2 if the user EXPLICITLY asks for bioisosteres or properties!

When user EXPLICITLY requests bioisosteres or properties for specific groups/sections:

STEP 1: Identify which groups the user wants to query
  - User may say: "find bioisosteres for the benzene ring"
  - Or: "find bioisosteres for group A and group B" (query EACH separately)
  - Or: "check all the rings" (query each ring from annotation)

STEP 2: Extract SMILES for EACH target group
  - Look at the previously annotated molecule's matches
  - For each referenced group, extract its SMILES representation
  - If user references by name (e.g., "benzene"), find matching annotation
  - If user references by index, use that match's atom_indices

STEP 3: Emit query set (REQUIRED before calling bioisostere/ADMET tools)

{"type":"query_set_structured","query_set":{"original_molecule_smiles":"<PARENT_SMILES>","targets":[{"smiles":"c1ccccc1","source_group":"benzene ring","atom_indices":[3,4,5,6,7,8]},{"smiles":"c1ccncc1","source_group":"pyridine ring","atom_indices":[10,11,12,13,14,15]}],"query_type":"bioisostere","reasoning":"User requested bioisosteres for the benzene and pyridine rings"}}

  - query_type: "bioisostere", "admet", or "combined"
  - targets: List of ALL groups to query (one entry per group)
  - Each target has: smiles, source_group name, and atom_indices from annotation

═══════════════════════════════════════════════════════════════════════════════
PHASE 3: BIOISOSTERE SCANNING (ONLY when user explicitly asks for bioisosteres)
═══════════════════════════════════════════════════════════════════════════════

⚠️ ONLY proceed to Phase 3 if the user EXPLICITLY asks for bioisosteres!
   Examples: "find bioisosteres", "suggest replacements", "scan for alternatives"
   Do NOT call this automatically after molecule annotation!

IMPORTANT: You must call scan_for_bioisosteres SEPARATELY for EACH target in the query set.
If user asked for 3 groups, make 3 separate MCP tool calls.

For EACH target in the query set:

STEP 1: Call the bioisostere scanner for this target
  - Call scan_for_bioisosteres(query_smiles=target.smiles, top_k=20, min_similarity=0.3)

STEP 2: Summarize results conversationally (NO JSON NEEDED)
  - The system automatically extracts results from MCP tool calls
  - DO NOT emit raw JSON - this wastes tokens!
  - Summarize findings naturally, e.g.:
    "I found 15 bioisosteres for the benzene ring. The top candidates include pyridine 
    (similarity 0.85), thiophene (0.72), and pyrrole (0.68). Would you like me to 
    get ADMET properties for any of these?"

STEP 3: Repeat for each target
  - If user asked for multiple groups, make separate tool calls for each
  - Summarize each result set conversationally

═══════════════════════════════════════════════════════════════════════════════
PHASE 4: ADMET PROPERTY PREDICTION (ONLY when user explicitly asks for properties)
═══════════════════════════════════════════════════════════════════════════════

⚠️ ONLY proceed to Phase 4 if the user EXPLICITLY asks for ADMET/properties!
   Examples: "predict properties", "check ADMET", "get solubility", "check drug-likeness"
   Do NOT call this automatically after molecule annotation!

When user EXPLICITLY requests properties (ADMET/pharmacokinetics):

STEP 1: Call the ADMET predictor for EACH SMILES
  - Call predict_admet_properties(smiles) separately for each molecule/fragment

STEP 2: Summarize ADMET results conversationally (NO JSON NEEDED)
  - The system automatically extracts results from MCP tool calls
  - DO NOT emit raw JSON with all predictions - this wastes tokens!
  - Summarize KEY findings naturally, e.g.:
    "The ADMET predictions look promising: good intestinal absorption (HIA: 0.95), 
    acceptable solubility, but watch out for moderate hERG liability (0.35). 
    The compound shows low CYP inhibition risk."
  - Focus on drug-likeness implications, not raw numbers
  - The frontend receives full prediction data automatically

═══════════════════════════════════════════════════════════════════════════════
CONVERSATIONAL GUIDELINES - ONLY CALL TOOLS WHEN EXPLICITLY NEEDED
═══════════════════════════════════════════════════════════════════════════════

WHEN YOU MUST USE TOOLS (MCP calls required):
  - User mentions ANY molecule by name (aspirin, caffeine, ibuprofen, etc.) → ONLY call get_smiles_from_name + annotate_molecule, then STOP
  - User asks for SMILES of anything → ONLY call get_smiles_from_name + annotate_molecule, then STOP
  - User provides a compound identifier → ONLY call convert_identifier_to_smiles + annotate_molecule, then STOP
  - User asks about structure/functional groups → ONLY call annotate_molecule, then STOP
  - User EXPLICITLY wants bioisosteres → call scan_for_bioisosteres (but ONLY if they ask!)
  - User EXPLICITLY wants ADMET/properties → call predict_admet_properties (but ONLY if they ask!)
  - User EXPLICITLY wants patent check → call check_molecule_patent (but ONLY if they ask!)

⚠️ CRITICAL: Do NOT proactively call tools the user didn't ask for!
   - If user asks for "aspirin", do NOT automatically check patents
   - If user asks for SMILES, do NOT automatically run ADMET
   - If user asks for structure, do NOT automatically scan bioisosteres
   - ONLY call the exact tools needed for their specific request

WHEN TO JUST TALK (no MCP call needed):
  - User asks a GENERAL chemistry question (no specific molecule mentioned)
  - User asks about data ALREADY returned by tools in this conversation
  - User asks for clarification on previous results
  - User asks which bioisostere to select from results already shown
  - User asks about mechanisms, drug design principles, SAR concepts
  - User wants to compare options already shown

WHEN TO GUIDE THE USER:
  - If user has bioisostere results but hasn't selected one:
    "I found several bioisosteres for the benzene ring. Would you like me to:
     - Get ADMET properties for any of these candidates?
     - Help you select based on specific criteria (e.g., lower logP, better solubility)?
     - Scan another ring in the molecule?"
  
  - If user's request is unclear:
    "I see several rings in this molecule. Which would you like to explore?
     1. The benzene ring (atoms 3-8)
     2. The pyridine ring (atoms 10-15)
     Or I can scan all of them."

STAY IN SCOPE - DRUG DEVELOPMENT ONLY:
  - You are a biochemist copilot for drug discovery and molecule design
  - Politely decline off-topic requests:
    "I'm specialized in drug development and molecule design. I can help with:
     - Analyzing molecular structures
     - Finding bioisosteric replacements
     - Predicting ADMET properties
     - Discussing SAR and medicinal chemistry
     Is there something in these areas I can help with?"
  
  - Do NOT discuss: politics, general coding, recipes, entertainment, etc.
  - Redirect firmly but kindly back to chemistry/pharma topics

SCIENTIFIC DISCUSSION IS WELCOME:
  - Explain WHY certain bioisosteres work (shape, electronics, H-bonding)
  - Discuss SAR implications of ring replacements
  - Explain ADMET property tradeoffs (e.g., "higher logP may improve permeability but hurt solubility")
  - Help interpret results in medicinal chemistry context
  - These discussions do NOT require tool calls

═══════════════════════════════════════════════════════════════════════════════
CRITICAL RULES
═══════════════════════════════════════════════════════════════════════════════

1. NEVER FABRICATE SMILES OR TOOL DATA
   - If user asks for SMILES → CALL get_smiles_from_name() - NEVER type SMILES yourself
   - If user mentions a molecule name → CALL get_smiles_from_name() first
   - After getting SMILES → ALWAYS call annotate_molecule() to provide full analysis
   - All bioisosteres must come from scan_for_bioisosteres (ONLY when user asks!)
   - All properties must come from predict_admet_properties (ONLY when user asks!)
   - Patent info must come from check_molecule_patent (ONLY when user asks!)
   - You CAN discuss data already returned by tools without re-calling

2. MOLECULE NAME = EXACTLY 2 TOOL CALLS (NO MORE, NO LESS)
   - "aspirin" → get_smiles_from_name("aspirin") + annotate_molecule(result) → STOP
   - "give me SMILES for X" → get_smiles_from_name("X") + annotate_molecule(result) → STOP
   - "what is the structure of Y" → get_smiles_from_name("Y") + annotate_molecule(result) → STOP
   - NEVER just respond with text containing SMILES without tool calls
   - NEVER call additional tools (patent, bioisostere, ADMET) unless user explicitly asks

3. DO NOT EMIT RAW JSON - THE SYSTEM HANDLES IT
   - The system automatically extracts structured data from MCP tool results
   - DO NOT regurgitate tool outputs as JSON - this wastes tokens!
   - Instead, summarize findings in natural language
   - The frontend receives full structured data via MCP, not your text output

4. TOOL ORDER MATTERS
   - name → SMILES → annotation → query set → bioisosteres/ADMET
   - Never skip the annotation step - always annotate after getting SMILES
   - annotate_molecule ONLY accepts SMILES (not names, not IDs)

5. ONE MCP CALL PER TARGET
   - If user wants bioisosteres for 3 groups → make 3 scan_for_bioisosteres calls
   - If user wants ADMET for 2 molecules → make 2 predict_admet_properties calls
   - Emit one structured JSON result after EACH tool call

6. BE A HELPFUL COPILOT
   - Guide users through the workflow naturally
   - Suggest next steps when appropriate
   - Answer scientific questions conversationally
   - Don't be robotic - be a knowledgeable colleague

═══════════════════════════════════════════════════════════════════════════════
AVAILABLE TOOLS
═══════════════════════════════════════════════════════════════════════════════

Phase 1 MCP (drugdiscovery_mcp):
  - get_smiles_from_name(chemical_name: str) → {"smiles": "..."}
  - convert_identifier_to_smiles(identifier: str) → {"smiles_list": [...], "count": N}
  - annotate_molecule(smiles: str) → {"name", "svg", "smiles", "matches": [...]}
  - scan_for_bioisosteres(query_smiles, top_k, min_similarity, ...) → {"results": [...]}
  - check_molecule_patent(smiles: str) → {"is_patented", "patents": [...]}

Phase 2 MCP (admet_mcp):
  - predict_admet_properties(smiles: str) → {"smiles", "predictions": {...}}
"""


class ChatAgent:
    """
    Chat agent using OpenAI's native MCP tool support.
    Connects to running MCP server via SSE (https://f38a55f3121b.ngrok-free.app by default).
    
    Prerequisites:
    - Run the MCP server: python backend/mcp_start.py
    - Set OPENAI_API_KEY environment variable
    """

    def __init__(self, model: str | None = None, mcp_server_url: str | None = None):
        self.model = model or os.getenv("DEFAULT_GPT_MODEL", DEFAULT_MODEL)
        self.client = OpenAI(api_key=os.getenv("OPENAI_API_KEY"))
        self.mcp_server_url_phase1 = _get_mcp_server_url("phase1_pipeline_mcp")
        self.mcp_server_url_phase2 = _get_mcp_server_url("phase2_pipeline_mcp")

        _dbg_print(_COLOR_GREEN, f"[AGENT] Initialized (model={self.model}, mcp_url_phase1={self.mcp_server_url_phase1}, mcp_url_phase2={self.mcp_server_url_phase2})")
    
    def chat(self, user_message: str, prior: Optional[List[Dict[str, str]]] = None) -> Dict[str, Any]:
        """
        Non-streaming chat using OpenAI's responses.create() with MCP tools.
        
        Args:
            user_message: User's message
            prior: Optional conversation history
        
        Returns:
            Dict containing:
            - text: Assistant's response text
            - tool_calls: List of tool calls made (if any)
        """
        _dbg_print(_COLOR_BLUE, f"[AGENT] chat request: {user_message[:100]}")
        
        # Build full conversation
        messages = [{"role": "system", "content": CHEM_DESIGN_SYSTEM_PROMPT}]
        if prior:
            messages.extend(prior)
        messages.append({"role": "user", "content": user_message})
        
        # Call OpenAI with MCP server
        _dbg_print(_COLOR_BLUE, f"[AGENT] Calling OpenAI with MCP server: {self.mcp_server_url_phase1} and {self.mcp_server_url_phase2}")
        
        try:
            resp = self.client.responses.create(
                model=self.model,
                tools=[
                    {
                        "type": "mcp",
                        "server_label": "phase1_annotation_pipeline_mcp",
                        "server_description": "Drug discovery tools: molecule annotation, identifier conversion, bioisostere scanning",
                        "server_url": self.mcp_server_url_phase1,
                        "require_approval": "never",
                    },
                    {
                        "type": "mcp",
                        "server_label": "phase2_pipeline_mcp",
                        "server_description": "Drug discovery tools: property prediction for new combined SMILES",
                        "server_url": self.mcp_server_url_phase2,
                        "require_approval": "never",
                    },
                ],
                input="\n".join(f"{m['role']}: {m['content']}" for m in messages),
            )
            
            # Extract tool calls from the response
            # OpenAI's Responses API with MCP stores results in resp.output items
            tool_calls = []
            
            # Debug: Log what attributes the response has
            resp_attrs = [a for a in dir(resp) if not a.startswith('_')]
            _dbg_print(_COLOR_BLUE, f"[AGENT] Response attributes: {resp_attrs}")
            
            # Deep debug: Try to serialize the response to see its structure
            try:
                if hasattr(resp, 'model_dump'):
                    resp_dict = resp.model_dump()
                    import json
                    _dbg_print(_COLOR_BLUE, f"[AGENT] Full response keys: {list(resp_dict.keys())}")
                    # Log output items if present
                    if 'output' in resp_dict:
                        _dbg_print(_COLOR_BLUE, f"[AGENT] Output has {len(resp_dict['output'])} items")
                        for i, item in enumerate(resp_dict['output'][:5]):  # First 5 items
                            item_type = item.get('type', 'unknown')
                            item_keys = list(item.keys()) if isinstance(item, dict) else []
                            _dbg_print(_COLOR_BLUE, f"[AGENT]   [{i}] type={item_type}, keys={item_keys}")
            except Exception as e:
                _dbg_print(_COLOR_RED, f"[AGENT] Could not serialize response: {e}")
            
            # Try to extract from resp.output (Responses API structure)
            if hasattr(resp, 'output') and resp.output:
                _dbg_print(_COLOR_BLUE, f"[AGENT] Found resp.output with {len(resp.output)} items")
                for idx, item in enumerate(resp.output):
                    item_type = getattr(item, 'type', None) or (item.get('type') if isinstance(item, dict) else None)
                    _dbg_print(_COLOR_BLUE, f"[AGENT]   Item {idx+1}: type={item_type}")
                    
                    # Look for MCP tool call results
                    if item_type in ('mcp_call', 'function_call', 'tool_use', 'mcp_tool_use'):
                        tool_name = getattr(item, 'name', None) or (item.get('name') if isinstance(item, dict) else 'unknown')
                        tool_args = getattr(item, 'arguments', {}) or (item.get('arguments', {}) if isinstance(item, dict) else {})
                        tool_result = getattr(item, 'output', None) or getattr(item, 'result', None) or (item.get('output') or item.get('result') if isinstance(item, dict) else None)
                        
                        _dbg_print(_COLOR_GREEN, f"[AGENT]   Found MCP call: {tool_name}")
                        
                        # Parse result if it's a string
                        if isinstance(tool_result, str):
                            try:
                                import json
                                tool_result = json.loads(tool_result)
                            except:
                                pass
                        
                        tool_calls.append({
                            "tool_name": tool_name,
                            "tool_args": tool_args if isinstance(tool_args, dict) else {},
                            "result": tool_result if isinstance(tool_result, dict) else {}
                        })
            
            # Fallback: Try resp.tool_calls (older API structure)
            if not tool_calls and hasattr(resp, 'tool_calls') and resp.tool_calls:
                _dbg_print(_COLOR_BLUE, f"[AGENT] Using resp.tool_calls fallback: {len(resp.tool_calls)} calls")
                for idx, tc in enumerate(resp.tool_calls):
                    tool_name = tc.get('name', 'unknown') if isinstance(tc, dict) else getattr(tc, 'name', 'unknown')
                    _dbg_print(_COLOR_BLUE, f"[AGENT]   Tool {idx+1}: {tool_name}")
                    
                    tool_call_info = {
                        "tool_name": tool_name,
                        "tool_args": tc.get('arguments', {}) if isinstance(tc, dict) else getattr(tc, 'arguments', {}),
                        "result": tc.get('result', {}) if isinstance(tc, dict) else getattr(tc, 'result', {})
                    }
                    tool_calls.append(tool_call_info)
            
            if tool_calls:
                _dbg_print(_COLOR_GREEN, f"[AGENT] ✅ Extracted {len(tool_calls)} tool call(s)")
            else:
                _dbg_print(_COLOR_RED, f"[AGENT] ⚠️  NO TOOL CALLS DETECTED - Model may be fabricating response!")
            
            _dbg_print(_COLOR_GREEN, f"[AGENT] response: {resp.output_text[:200]}...")
            
            return {
                "text": resp.output_text,
                "tool_calls": tool_calls if tool_calls else None
            }
            
        except AttributeError as e:
            _dbg_print(_COLOR_RED, f"[AGENT] ⚠️  API Error: {e}")
            _dbg_print(_COLOR_RED, f"[AGENT] responses.create() may not be available - check OpenAI API version")
            raise
    
    def guide(self, user_message: str, prior: Optional[List[Dict[str, str]]] = None) -> Dict[str, Any]:
        """Backward-compatible helper: molecule workflow guide."""
        return self.chat(user_message, prior)
    
    def copilot(self, user_message: str, prior: Optional[List[Dict[str, str]]] = None) -> Dict[str, Any]:
        """Backward-compatible helper: biochemist copilot."""
        return self.chat(user_message, prior)
