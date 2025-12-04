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
    # This is a more robust approach that handles nested braces
    json_pattern = rf'\{{\s*["\']type["\']\s*:\s*["\'](?:{type_pattern})["\']'
    
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
    
    # Handle multiple bioisostere or admet results
    if len(structured_objects) > 1:
        types = [obj.get("type") for obj in structured_objects]
        if all(t == "bioisostere_structured" for t in types):
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
    - multi_bioisostere: Multiple bioisostere results (e.g., for "A OR B")
    - multi_admet: Multiple ADMET predictions
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
        
        # Extract any structured data embedded in the response
        cleaned_text, structured_objects, response_type = extract_structured_data(response_text)
        
        _debug_print(f"Response type: {response_type}", _COLOR_GREEN)
        _debug_print(f"Structured objects found: {len(structured_objects)}", _COLOR_GREEN)
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
        _debug_print(f"bioisostere_data: {'YES' if response.bioisostere_data else 'NO'}", _COLOR_RED)
        _debug_print(f"admet_data: {'YES' if response.admet_data else 'NO'}", _COLOR_RED)
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
    Populate the typed data fields based on structured objects and response type.
    """
    try:
        if response_type == "annotation" and structured_objects:
            obj = structured_objects[0]
            if "annotation" in obj:
                response.molecule_data = MoleculeAnnotation(**obj["annotation"])
            elif "smiles" in obj and "matches" in obj:
                response.molecule_data = MoleculeAnnotation(**obj)
                
        elif response_type == "query_set" and structured_objects:
            obj = structured_objects[0]
            if "query_set" in obj:
                response.query_set_data = SmilesQuerySet(**obj["query_set"])
            else:
                response.query_set_data = SmilesQuerySet(**obj)
                
        elif response_type == "bioisosteres" and structured_objects:
            obj = structured_objects[0]
            response.bioisostere_data = BioisostereScanResult(
                query_smiles=obj.get("query_smiles", ""),
                num_results=len(obj.get("results", [])),
                results=obj.get("results", [])
            )
            
        elif response_type == "admet" and structured_objects:
            obj = structured_objects[0]
            response.admet_data = ADMETPropertyResult(
                smiles=obj.get("smiles", ""),
                predictions=obj.get("predictions", {})
            )
            
        elif response_type == "multi_bioisostere" and structured_objects:
            # TODO: Build MultiQueryBioisostereResult from multiple objects
            pass
            
        elif response_type == "multi_admet" and structured_objects:
            # TODO: Build MultiQueryADMETResult from multiple objects
            pass
            
    except Exception as e:
        # Log but don't fail - structured_data dict is still available
        import logging
        logging.warning(f"Failed to populate typed fields: {e}")
    
    return response