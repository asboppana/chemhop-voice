"""
Voice backend with MCP command queueing system.

Commands the LLM can perform:
- Create a new protein
- Give you a list of similar proteins

In the model's system prompt it is informed of:
- The conversation history
- What commands it is able to and not able to perform
- It must continue the conversation even when commands are running or queued
"""

import logging
from elevenlabs.client import ElevenLabs
from src.config.settings import get_settings
from src.services.elevenlabs.knowledge_base import get_knowledge_base
from src.services.prompts import VOICE_AGENT_SYSTEM_PROMPT, VOICE_AGENT_FIRST_MESSAGE


def create_elevenlabs_agent(name: str) -> str:
    settings = get_settings()
    knowledge_base = get_knowledge_base()
    elevenlabs_client = ElevenLabs(api_key=settings.elevenlabs_api_key)
    response = elevenlabs_client.conversational_ai.agents.create(
        name=name,
        conversation_config={
            "turn": {
                "soft_timeout_config": {"timeout_seconds": 3.0, "message": "Hmmm, well that's interesting..."}
            },
            "tts": {
                "voice_id": "21m00Tcm4TlvDq8ikWAM",
                "model_id": "eleven_turbo_v2",
                "speed": 1.0,
                "stability": 0.8
            },
            "conversation": {
                "max_duration_seconds": 600,
                "client_events": [
                    "agent_response", 
                    "agent_response_correction",
                    "agent_response_metadata",
                    "audio",
                    "interruption",
                    "vad_score",
                    "user_transcript",
                ]
            },
            "agent": {
                "first_message": VOICE_AGENT_FIRST_MESSAGE,
                "prompt": {
                    "prompt": VOICE_AGENT_SYSTEM_PROMPT,
                    "llm": "claude-sonnet-4-5",
                    "temperature": 0.7,
                    "knowledge_base": knowledge_base
                }
            }
        }
    )

    logging.info(f"Agent created with ID: {response.agent_id}")
    return response.agent_id

if __name__ == "__main__":
    agent_id = create_elevenlabs_agent("drugdiscoveryagent")
