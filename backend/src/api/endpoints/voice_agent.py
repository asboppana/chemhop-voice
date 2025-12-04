"""
Voice Agent endpoints.
Provides access to ElevenLabs voice agent configuration and chat functionality.
"""
import json
import re
from typing import List, Dict, Any, Optional, Tuple
from fastapi import APIRouter, HTTPException

from src.api.models.chat import ChatRequest
from src.api.models.mcp_models import (
    AgentChatResponse,
    MCPToolCall,
    MoleculeAnnotation,
    SmilesQuerySet,
    BioisostereScanResult,
    ADMETPropertyResult,
    MultiQueryBioisostereResult,
    MultiQueryADMETResult,
)
from src.services.chatAgent.agent import ChatAgent

# Color codes for debug printing
_COLOR_RED = "\033[91m"
_COLOR_GREEN = "\033[92m"
_COLOR_YELLOW = "\033[93m"
_COLOR_BLUE = "\033[94m"
_COLOR_RESET = "\033[0m"


def _debug_print(message: str, color: str = _COLOR_RED) -> None:
    """Print debug message in color."""
    print(f"{color}[VOICE_AGENT] {message}{_COLOR_RESET}", flush=True)


def _extract_from_mcp_results(tool_calls_raw: Optional[List[Dict[str, Any]]]) -> Optional[Dict[str, Any]]:
    """
    Extract structured data directly from MCP tool call results.
    
    This is the PREFERRED method - it uses Pydantic validation on actual MCP results
    instead of parsing JSON from the model's text output.
    
    Benefits:
    - No fragile regex/JSON parsing
    - Proper Pydantic validation with clear error messages
    - Uses actual MCP data, not model's regurgitated version
    - Saves tokens (model doesn't need to regenerate huge SVG strings)
    
    Returns:
        Dict with 'objects' (list of structured dicts) and 'response_type' (str),
        or None if no relevant tool calls found.
    """
    if not tool_calls_raw:
        return None
    
    structured_objects = []
    response_type = "text_only"
    
    for tc in tool_calls_raw:
        tool_name = tc.get("tool_name", "")
        result = tc.get("result", {})
        
        if not result or isinstance(result, str):
            continue
        
        # Handle annotate_molecule results
        if tool_name == "annotate_molecule":
            try:
                # Validate with Pydantic - this will raise clear errors if invalid
                annotation = MoleculeAnnotation(**result)
                
                # Build structured object in the expected format
                structured_obj = {
                    "type": "molecule_structured",
                    "smiles": annotation.smiles,
                    "annotation": annotation.model_dump()
                }
                structured_objects.append(structured_obj)
                
                if response_type == "text_only":
                    response_type = "annotation"
                    
                _debug_print(f"  ✅ annotate_molecule: validated with Pydantic ({len(annotation.matches)} matches)", _COLOR_GREEN)
                
            except Exception as e:
                _debug_print(f"  ❌ annotate_molecule validation failed: {e}", _COLOR_RED)
                # Fall through to let text parsing try
                continue
        
        # Handle scan_for_bioisosteres results
        elif tool_name == "scan_for_bioisosteres":
            try:
                # Validate with Pydantic
                scan_result = BioisostereScanResult(**result)
                
                structured_obj = {
                    "type": "bioisostere_structured",
                    "query_smiles": scan_result.query_smiles,
                    "source_group": scan_result.source_group,
                    "atom_indices": scan_result.atom_indices,
                    "results": [r.model_dump() for r in scan_result.results]
                }
                structured_objects.append(structured_obj)
                
                if response_type == "text_only":
                    response_type = "bioisosteres"
                elif response_type == "bioisosteres":
                    response_type = "multi_bioisostere"
                    
                _debug_print(f"  ✅ scan_for_bioisosteres: validated ({scan_result.num_results} results)", _COLOR_GREEN)
                
            except Exception as e:
                _debug_print(f"  ❌ scan_for_bioisosteres validation failed: {e}", _COLOR_RED)
                continue
        
        # Handle predict_admet_properties results
        elif tool_name == "predict_admet_properties":
            try:
                # Validate with Pydantic
                admet_result = ADMETPropertyResult(**result)
                
                structured_obj = {
                    "type": "admet_structured",
                    "smiles": admet_result.smiles,
                    "predictions": admet_result.predictions
                }
                structured_objects.append(structured_obj)
                
                if response_type == "text_only":
                    response_type = "admet"
                elif response_type == "admet":
                    response_type = "multi_admet"
                    
                _debug_print(f"  ✅ predict_admet_properties: validated ({len(admet_result.predictions)} properties)", _COLOR_GREEN)
                
            except Exception as e:
                _debug_print(f"  ❌ predict_admet_properties validation failed: {e}", _COLOR_RED)
                continue
    
    # Determine combined response types
    if len(structured_objects) > 1:
        types = [obj.get("type") for obj in structured_objects]
        unique_types = set(types)
        
        if "molecule_structured" in unique_types and "admet_structured" in unique_types:
            response_type = "molecule_with_admet"
        elif "molecule_structured" in unique_types and "bioisostere_structured" in unique_types:
            response_type = "molecule_with_bioisostere"
        elif len(unique_types) > 1:
            response_type = "combined"
    
    if not structured_objects:
        return None
    
    return {
        "objects": structured_objects,
        "response_type": response_type
    }


def _clean_json_from_text(text: str) -> str:
    """
    Remove any JSON objects from the text that match our structured types.
    
    When we extract data directly from MCP results, the model might still
    emit JSON in its response. We clean it out to avoid confusion.
    """
    type_pattern = "|".join(STRUCTURED_TYPES)
    json_pattern = rf'\{{"type"\s*:\s*"(?:{type_pattern})"'
    
    cleaned = text
    for match in re.finditer(json_pattern, text):
        start_pos = match.start()
        
        # Find the complete JSON object
        brace_count = 0
        end_pos = start_pos
        in_string = False
        escape_next = False
        
        for i, char in enumerate(text[start_pos:], start=start_pos):
            if escape_next:
                escape_next = False
                continue
            if char == '\\':
                escape_next = True
                continue
            if char == '"' and not escape_next:
                in_string = not in_string
                continue
            if in_string:
                continue
            if char == '{':
                brace_count += 1
            elif char == '}':
                brace_count -= 1
                if brace_count == 0:
                    end_pos = i + 1
                    break
        
        if end_pos > start_pos:
            json_str = text[start_pos:end_pos]
            cleaned = cleaned.replace(json_str, '').strip()
    
    # Clean up extra whitespace
    cleaned = re.sub(r'\n\s*\n', '\n\n', cleaned).strip()
    return cleaned


router = APIRouter(tags=["voice_agent"])

# Structured message type identifiers
STRUCTURED_TYPES = [
    "molecule_structured",
    "query_set_structured", 
    "bioisostere_structured",
    "admet_structured"
]


def extract_structured_data(response_text: str) -> Tuple[str, List[Dict[str, Any]], Optional[str]]:
    """
    Extract structured JSON messages from agent response text.
    
    The agent may embed JSON objects in its response like:
    {"type":"molecule_structured","smiles":"...","annotation":{...}}
    {"type":"query_set_structured","query_set":{...}}
    {"type":"bioisostere_structured","query_smiles":"...","results":[...]}
    {"type":"admet_structured","smiles":"...","predictions":{...}}
    
    This function:
    1. Finds and parses these JSON objects
    2. Removes them from the text
    3. Returns cleaned text, list of structured objects, and detected response type
    
    Args:
        response_text: Raw response from the agent
        
    Returns:
        tuple: (cleaned_text, list_of_structured_objects, response_type)
    """
    structured_objects = []
    cleaned_text = response_text
    response_type = "text_only"
    
    # Build pattern for all structured types
    type_pattern = "|".join(STRUCTURED_TYPES)
    
    # Pattern to match JSON objects with nested structures  
    # Matches: {"type":"molecule_structured"...
    json_pattern = rf'{{"type"\s*:\s*"(?:{type_pattern})"'
    
    # Find all starting positions of structured messages
    for match in re.finditer(json_pattern, response_text):
        start_pos = match.start()
        
        # Extract the complete JSON object by counting braces
        brace_count = 0
        end_pos = start_pos
        in_string = False
        escape_next = False
        
        for i, char in enumerate(response_text[start_pos:], start=start_pos):
            if escape_next:
                escape_next = False
                continue
            if char == '\\':
                escape_next = True
                continue
            if char == '"' and not escape_next:
                in_string = not in_string
                continue
            if in_string:
                continue
            if char == '{':
                brace_count += 1
            elif char == '}':
                brace_count -= 1
                if brace_count == 0:
                    end_pos = i + 1
                    break
        
        if end_pos > start_pos:
            json_str = response_text[start_pos:end_pos]
            try:
                json_obj = json.loads(json_str)
                structured_objects.append(json_obj)
                # Update response type based on first structured message
                if response_type == "text_only":
                    msg_type = json_obj.get("type", "")
                    if msg_type == "molecule_structured":
                        response_type = "annotation"
                    elif msg_type == "query_set_structured":
                        response_type = "query_set"
                    elif msg_type == "bioisostere_structured":
                        response_type = "bioisosteres"
                    elif msg_type == "admet_structured":
                        response_type = "admet"
                # Remove from text
                cleaned_text = cleaned_text.replace(json_str, '').strip()
            except json.JSONDecodeError:
                continue
    
    # Handle multiple results
    if len(structured_objects) > 1:
        types = [obj.get("type") for obj in structured_objects]
        unique_types = set(types)
        
        # Check for mixed types (e.g., molecule + admet)
        if len(unique_types) > 1:
            # Determine combined response type
            if "molecule_structured" in unique_types and "admet_structured" in unique_types:
                response_type = "molecule_with_admet"
            elif "molecule_structured" in unique_types and "bioisostere_structured" in unique_types:
                response_type = "molecule_with_bioisostere"
            else:
                response_type = "combined"
        # All same type
        elif all(t == "bioisostere_structured" for t in types):
            response_type = "multi_bioisostere"
        elif all(t == "admet_structured" for t in types):
            response_type = "multi_admet"
    
    # Clean up extra whitespace
    cleaned_text = re.sub(r'\n\s*\n', '\n\n', cleaned_text).strip()
    
    return cleaned_text, structured_objects, response_type


@router.post("/chat", response_model=AgentChatResponse)
async def chat_with_agent(request: ChatRequest) -> AgentChatResponse:
    """
    Chat with the biochemist copilot agent.
    
    This endpoint processes chat messages using the ChatAgent which integrates
    with the MCP servers for drug discovery tools:
    
    Phase 1 MCP (drugdiscovery_mcp):
    - annotate_molecule: SMILES → MoleculeAnnotation
    - get_smiles_from_name: chemical_name → SMILES
    - convert_identifier_to_smiles: identifier → SMILES
    - check_molecule_patent: SMILES → PatentCheckResult
    - scan_for_bioisosteres: SMILES → BioisostereScanResult
    
    Phase 2 MCP (admet_mcp):
    - predict_admet_properties: SMILES → ADMETPropertyResult
    
    Response Types:
    - annotation: Molecule annotation with functional groups
    - query_set: Set of SMILES extracted from user intent
    - bioisosteres: Single bioisostere scan result
    - admet: Single ADMET prediction
    - molecule_with_admet: Molecule annotation + ADMET prediction (both fields populated)
    - molecule_with_bioisostere: Molecule annotation + bioisostere results (both fields populated)
    - multi_bioisostere: Multiple bioisostere results (e.g., for "A OR B")
    - multi_admet: Multiple ADMET predictions
    - combined: Other combinations of structured data
    - text_only: Pure text response (no structured data)
    
    Args:
        request: ChatRequest with messages list and optional model
    
    Returns:
        AgentChatResponse: Structured response with text and typed data fields
    """
    try:
        # Extract conversation history
        messages = request.messages
        if not messages:
            raise HTTPException(
                status_code=400,
                detail="Messages list cannot be empty"
            )
        
        # Last message should be from user
        user_message = messages[-1].get("content", "")
        if not user_message:
            raise HTTPException(
                status_code=400,
                detail="Last message must have content"
            )
        
        # Prior messages are everything except the last one
        prior = messages[:-1] if len(messages) > 1 else None
        
        # Initialize agent with optional model
        _debug_print(f"Initializing agent with model: {request.model}", _COLOR_BLUE)
        agent = ChatAgent(model=request.model)
        
        # Get response from agent (synchronous call)
        # Returns dict with "text" and "tool_calls" fields
        _debug_print(f"Sending message to agent: {user_message[:100]}...", _COLOR_BLUE)
        agent_response = agent.chat(user_message=user_message, prior=prior)
        response_text = agent_response.get("text", "")
        tool_calls_raw = agent_response.get("tool_calls")
        
        _debug_print(f"Raw response text ({len(response_text)} chars): {response_text[:500]}...", _COLOR_YELLOW)
        _debug_print(f"Tool calls: {len(tool_calls_raw) if tool_calls_raw else 0}", _COLOR_YELLOW)
        
        # OPTIMIZATION: Extract structured data directly from MCP tool results
        # This avoids parsing JSON from model output and uses Pydantic validation
        mcp_structured_data = _extract_from_mcp_results(tool_calls_raw)
        
        if mcp_structured_data:
            _debug_print(f"✅ Extracted structured data directly from MCP results!", _COLOR_GREEN)
            # Use MCP results directly - no need to parse model's text
            structured_objects = mcp_structured_data["objects"]
            response_type = mcp_structured_data["response_type"]
            # Clean the response text of any JSON the model might have emitted
            cleaned_text = _clean_json_from_text(response_text)
        else:
            # Fallback: try to extract from model text (for non-MCP responses)
            cleaned_text, structured_objects, response_type = extract_structured_data(response_text)
        
        _debug_print(f"Response type: {response_type}", _COLOR_GREEN)
        _debug_print(f"Structured objects found: {len(structured_objects)}", _COLOR_GREEN)
        for i, obj in enumerate(structured_objects):
            obj_type = obj.get("type", "unknown")
            _debug_print(f"  Object {i+1}: type={obj_type}", _COLOR_GREEN)
            if obj_type == "molecule_structured":
                has_annotation = "annotation" in obj
                has_matches = "matches" in obj.get("annotation", {}) if has_annotation else "matches" in obj
                num_matches = len(obj.get("annotation", {}).get("matches", [])) if has_annotation else len(obj.get("matches", []))
                _debug_print(f"    → has_annotation={has_annotation}, matches={num_matches}", _COLOR_GREEN)
        _debug_print(f"Cleaned text: {cleaned_text[:300]}..." if cleaned_text else "No cleaned text", _COLOR_GREEN)
        
        # Build structured response
        structured_data = None
        if structured_objects:
            structured_data = {"embedded_objects": structured_objects}
            for i, obj in enumerate(structured_objects):
                obj_type = obj.get("type", "unknown")
                _debug_print(f"  Object {i+1}: type={obj_type}", _COLOR_GREEN)
        
        # Convert tool calls to MCPToolCall format if present
        tool_calls = None
        if tool_calls_raw:
            tool_calls = [
                MCPToolCall(
                    tool_name=tc["tool_name"],
                    server=_infer_server_from_tool(tc["tool_name"]),
                    tool_args=tc["tool_args"],
                    result=tc["result"]
                )
                for tc in tool_calls_raw
            ]
        
        # Build response with typed fields based on response_type
        response = AgentChatResponse(
            response=cleaned_text if structured_objects else response_text,
            status="success",
            response_type=response_type,
            tool_calls=tool_calls,
            structured_data=structured_data,
        )
        
        # Populate typed fields based on structured objects
        if structured_objects:
            response = _populate_typed_fields(response, structured_objects, response_type)
        
        _debug_print(f"=== FINAL RESPONSE ===", _COLOR_RED)
        _debug_print(f"response_type: {response.response_type}", _COLOR_RED)
        _debug_print(f"response text: {response.response[:200]}..." if response.response else "None", _COLOR_RED)
        _debug_print(f"molecule_data: {'YES' if response.molecule_data else 'NO'}", _COLOR_RED)
        if response.molecule_data:
            num_matches = len(response.molecule_data.matches)
            matches_with_svg = sum(1 for m in response.molecule_data.matches if m.svg)
            _debug_print(f"  → {num_matches} matches, {matches_with_svg} with SVG ({'✅' if matches_with_svg > 0 else '❌'})", _COLOR_RED)
        _debug_print(f"query_set_data: {'YES' if response.query_set_data else 'NO'}", _COLOR_RED)
        _debug_print(f"bioisostere_data: {'YES' if response.bioisostere_data else 'NO'}", _COLOR_RED)
        _debug_print(f"admet_data: {'YES' if response.admet_data else 'NO'}", _COLOR_RED)
        _debug_print(f"multi_bioisostere_data: {'YES' if response.multi_bioisostere_data else 'NO'}", _COLOR_RED)
        _debug_print(f"multi_admet_data: {'YES' if response.multi_admet_data else 'NO'}", _COLOR_RED)
        _debug_print(f"======================", _COLOR_RED)
        
        return response
        
    except HTTPException:
        raise
    except Exception as e:
        import traceback
        traceback.print_exc()
        raise HTTPException(
            status_code=500,
            detail=f"Error processing chat request: {str(e)}"
        )


def _infer_server_from_tool(tool_name: str) -> str:
    """Infer which MCP server a tool belongs to."""
    phase2_tools = {"predict_admet_properties"}
    return "phase2" if tool_name in phase2_tools else "phase1"


def _populate_typed_fields(
    response: AgentChatResponse, 
    structured_objects: List[Dict[str, Any]], 
    response_type: str
) -> AgentChatResponse:
    """
    Populate the typed data fields based on structured objects.
    This now handles multiple different types in a single response.
    """
    try:
        # Process each structured object and populate the appropriate field
        # Use independent if statements (not elif) to handle multiple types
        for obj in structured_objects:
            obj_type = obj.get("type", "")
            
            if obj_type == "molecule_structured":
                # Populate molecule_data
                _debug_print(f"Processing molecule_structured object", _COLOR_GREEN)
                if "annotation" in obj:
                    # Merge top-level smiles into annotation if missing
                    annotation_data = obj["annotation"].copy()
                    if "smiles" not in annotation_data and "smiles" in obj:
                        annotation_data["smiles"] = obj["smiles"]
                    
                    # Debug: Check if matches have SVG
                    num_matches = len(annotation_data.get("matches", []))
                    matches_with_svg = sum(1 for m in annotation_data.get("matches", []) if m.get("svg"))
                    _debug_print(f"  Annotation: {num_matches} matches, {matches_with_svg} with SVG", _COLOR_GREEN)
                    
                    response.molecule_data = MoleculeAnnotation(**annotation_data)
                    _debug_print(f"  ✅ molecule_data populated successfully", _COLOR_GREEN)
                elif "smiles" in obj and "matches" in obj:
                    _debug_print(f"  Using flat structure (no annotation wrapper)", _COLOR_GREEN)
                    response.molecule_data = MoleculeAnnotation(**obj)
                    
            elif obj_type == "query_set_structured":
                # Populate query_set_data
                if "query_set" in obj:
                    response.query_set_data = SmilesQuerySet(**obj["query_set"])
                else:
                    response.query_set_data = SmilesQuerySet(**obj)
                    
            elif obj_type == "bioisostere_structured":
                # Populate bioisostere_data (single result)
                response.bioisostere_data = BioisostereScanResult(
                    query_smiles=obj.get("query_smiles", ""),
                    source_group=obj.get("source_group"),
                    atom_indices=obj.get("atom_indices"),
                    num_results=len(obj.get("results", [])),
                    results=obj.get("results", [])
                )
                
            elif obj_type == "admet_structured":
                # Populate admet_data (single result)
                response.admet_data = ADMETPropertyResult(
                    smiles=obj.get("smiles", ""),
                    predictions=obj.get("predictions", {})
                )
        
        # Handle multi-query results (all same type)
        if response_type == "multi_bioisostere":
            # TODO: Build MultiQueryBioisostereResult from multiple bioisostere objects
            pass
            
        elif response_type == "multi_admet":
            # TODO: Build MultiQueryADMETResult from multiple admet objects
            pass
            
    except Exception as e:
        # Log but don't fail - structured_data dict is still available
        import logging
        logging.warning(f"Failed to populate typed fields: {e}")
    
    return response