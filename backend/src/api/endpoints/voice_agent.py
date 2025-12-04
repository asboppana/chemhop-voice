"""
Voice Agent endpoints.
Provides access to ElevenLabs voice agent configuration.
"""
from fastapi import APIRouter, Request, HTTPException

from src.api.models.voice import VoiceAgentResponse

router = APIRouter(tags=["voice_agent"])


@router.get("/voice-agent", response_model=VoiceAgentResponse)
async def get_voice_agent_id(request: Request):
    """
    Get the ElevenLabs voice agent ID.
    
    This endpoint returns the agent ID that was created during application startup.
    The agent ID is used by the frontend to initialize the ElevenLabs conversation widget.
    """
    agent_id = getattr(request.app.state, "voice_agent_id", None)
    print(f"\033[91m{agent_id}\033[0m")
    if not agent_id:
        raise HTTPException(
            status_code=503,
            detail="Voice agent not initialized. Please try again later."
        )
    
    return VoiceAgentResponse(agent_id=agent_id)

