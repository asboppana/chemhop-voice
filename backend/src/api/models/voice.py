"""
Request and response models for voice agent endpoint.
"""
from pydantic import BaseModel

class VoiceAgentResponse(BaseModel):
    """Voice agent configuration response model."""
    agent_id: str
