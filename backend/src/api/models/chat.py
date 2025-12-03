"""
Request and response models for LLM chat endpoints.
"""
from typing import Any, Dict, List, Optional

from pydantic import BaseModel


class ChatRequest(BaseModel):
    """Payload for health chat.
    
    - messages: List of message dictionaries with role and content
    - model: Optional model name supported by LiteLLM (defaults from env or sane default)
    """
    messages: List[Dict[str, Any]]  # [{role: "user" | "assistant", content: str}]
    model: Optional[str] = None

