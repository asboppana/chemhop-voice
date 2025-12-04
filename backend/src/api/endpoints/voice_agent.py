"""
Voice Agent endpoints.
Provides access to ElevenLabs voice agent configuration and chat functionality.
"""
from fastapi import APIRouter, Request, HTTPException

from src.api.models.voice import VoiceAgentResponse
from src.api.models.chat import ChatRequest
from src.services.chatAgent.agent import ChatAgent

router = APIRouter(tags=["voice_agent"])


@router.post("/chat")
async def chat_with_agent(request: ChatRequest):
    """
    Chat with the biochemist copilot agent.
    
    This endpoint processes chat messages using the ChatAgent which integrates
    with the MCP server for drug discovery tools (molecule annotation, 
    bioisostere scanning, ADMET prediction, etc.).
    
    Args:
        request: ChatRequest with messages list and optional model
    
    Returns:
        dict: Response containing the agent's message
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
        agent = ChatAgent(model=request.model)
        
        # Get response from agent (synchronous call)
        response_text = agent.chat(user_message=user_message, prior=prior)
        
        return {
            "response": response_text,
            "status": "success"
        }
        
    except HTTPException:
        raise
    except Exception as e:
        import traceback
        traceback.print_exc()
        raise HTTPException(
            status_code=500,
            detail=f"Error processing chat request: {str(e)}"
        )