"""
Streaming LLM Agent with native MCP server integration.

This agent:
- Uses GPT-5 (or gpt-4o) with native MCP tool support
- Connects to our local MCP server (drugdiscovery_mcp) via stdio/SSE
- Exposes tools: convert_identifier_to_smiles, annotate_molecule, scan_for_bioisosteres, etc.
- System prompt guides scaffold hopping and bioisostere-driven molecule design
"""

from __future__ import annotations

import json
import logging
import os
from typing import Dict, List, Optional

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
    Set MCP_SERVER_URL env var to override.
    """
    # Use env var if set, otherwise default to local MCP SSE server
    # For ngrok: export MCP_SERVER_URL=https://your-ngrok-url.ngrok-free.app/sse
    if server_label == "phase1_pipeline_mcp":
        return "https://f38a55f3121b.ngrok-free.app/sse"
    elif server_label == "phase2_pipeline_mcp":
        return "https://kaycee-radiculose-novelistically.ngrok-free.dev/sse"
    else:   
        raise ValueError(f"Invalid server label: {server_label}")


CHEM_DESIGN_SYSTEM_PROMPT = (
    "You are a biochemist copilot guiding scaffold hopping and bioisostere-driven design.\n"
    "Operate in three phases with structured planning:\n"
    "  - ingest: Identify/normalize the molecule (identifier, which can be normal languge -> SMILES) and annotate functional groups\n"
    "  - design: Interpret user's modification intent (boolean expressions over annotated sections),\n"
    "            find ring bioisosteres for requested sections, and enumerate candidate replacements\n"
    "  - properties: Compute properties for candidate rings or provided modified molecules\n"
    "\n"
    "Rules:\n"
    "- Always prefer tools over guesses. When user provides a molecule identifier:\n"
    "  1. If it's a chemical name, ChEMBL ID (chembl:XXX), ChEBI ID (chebi:XXX), PubChem CID (pubchem:XXX),\n"
    "     or you're unsure if it's valid SMILES, call convert_identifier_to_smiles first\n"
    "  2. Then call annotate_molecule with the SMILES string (annotate_molecule only accepts SMILES format)\n"
    "  3. IMPORTANT: annotate_molecule requires a valid SMILES string - never pass identifiers directly\n"
    "\n"
    "- After you acquire both SMILES and annotation, emit a structured JSON message immediately:\n"
    "  Format: {\\\"type\\\":\\\"molecule_structured\\\",\\\"smiles\\\":\\\"<SMILES_STRING>\\\",\\\"annotation\\\":<FULL_ANNOTATION_OBJECT>}\n"
    "  The annotation object should contain: name, svg, matches (with atom_indices, trivial_name, etc.)\n"
    "\n"
    "- Users can specify modifications like: 'section A AND section B' or 'section A OR (section B AND section C)'.\n"
    "  Interpret 'sections' as specific annotated groups or ring substructures. If unclear, ask the user to select groups.\n"
    "- Boolean logic semantics:\n"
    "  AND = combine modifications across all referenced sections (Cartesian product).\n"
    "  OR = any of the referenced sections individually (set union of candidates).\n"
    "- For ring-focused replacements, extract a ring fragment SMILES from the annotated atom indices before scanning.\n"
    "  Use fragment_from_atoms if needed. Then call scan_for_bioisosteres for each target ring.\n"
    "- When the user asks for properties, call predict_admet_properties for each candidate SMILES provided\n"
    "  (ring or full molecule). If you only have ring candidates, note these are ring-level properties.\n"
    "- If no group/section is selected for modification, clearly ask the user to choose (e.g., by group index/name).\n"
    "- Be concise. Stream short guidance while tools run. Never fabricate tool outputs.\n"
)


# (No manual bridge or subprocess needed - OpenAI communicates with running MCP server via SSE)


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
    
    def chat(self, user_message: str, prior: Optional[List[Dict[str, str]]] = None) -> str:
        """
        Non-streaming chat using OpenAI's responses.create() with MCP tools.
        
        Args:
            user_message: User's message
            prior: Optional conversation history
        
        Returns:
            Assistant's response text
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
            
            # Debug: Log tool calls if present
            if hasattr(resp, 'tool_calls') and resp.tool_calls:
                _dbg_print(_COLOR_BLUE, f"[AGENT] Tool calls made: {len(resp.tool_calls)}")
                for idx, tc in enumerate(resp.tool_calls):
                    _dbg_print(_COLOR_BLUE, f"[AGENT]   Tool {idx+1}: {tc.get('name', 'unknown')}")
            else:
                _dbg_print(_COLOR_RED, f"[AGENT] ⚠️  NO TOOL CALLS DETECTED - Model may be fabricating response!")
            
            _dbg_print(_COLOR_GREEN, f"[AGENT] response: {resp.output_text[:200]}...")
            return resp.output_text
            
        except AttributeError as e:
            _dbg_print(_COLOR_RED, f"[AGENT] ⚠️  API Error: {e}")
            _dbg_print(_COLOR_RED, f"[AGENT] responses.create() may not be available - check OpenAI API version")
            raise
    
    def guide(self, user_message: str, prior: Optional[List[Dict[str, str]]] = None) -> str:
        """Backward-compatible helper: molecule workflow guide."""
        return self.chat(user_message, prior)
    
    def copilot(self, user_message: str, prior: Optional[List[Dict[str, str]]] = None) -> str:
        """Backward-compatible helper: biochemist copilot."""
        return self.chat(user_message, prior)
