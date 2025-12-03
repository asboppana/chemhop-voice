"""
LLM controller for health chat functionality.

Handles streaming LLM completions using OpenAI API.
"""
import logging
from typing import AsyncGenerator

from fastapi import HTTPException, status
from openai import AsyncOpenAI

from src.api.models.chat import ChatRequest
from src.config.settings import get_settings
from src.services.prompts import CHAT_SYSTEM_PROMPT, CHAT_END_PROMPT

logger = logging.getLogger(__name__)

# Default model to use when none is specified
DEFAULT_MODEL = "gpt-4-turbo"


class LLMController:
    """Controller for LLM chat operations."""

    def __init__(self):
        """Initialize LLM controller with OpenAI client."""
        settings = get_settings()
        self.client = AsyncOpenAI(api_key=settings.openai_api_key)

    def _validate_request(self, request: ChatRequest) -> None:
        """
        Validate health chat request.
        
        Args:
            request: ChatRequest with messages and optional model
            
        Raises:
            HTTPException 400: If messages are empty or invalid
        """
        if not request.messages:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail="messages cannot be empty"
            )
        
        # Get last user message content
        user_text = request.messages[-1].get("content", "")
        if not user_text:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail="health_info cannot be empty"
            )

    async def health_chat_stream(
        self, request: ChatRequest
    ) -> AsyncGenerator[str, None]:
        """
        Generate streaming health chat response.
        
        NOTE: Validation should be done BEFORE calling this method (in the endpoint)
        to ensure validation errors are raised before streaming starts.
        
        Args:
            request: ChatRequest with messages and optional model
            
        Yields:
            Text chunks as they are generated
        """
        # Set default model if none provided
        model = request.model or DEFAULT_MODEL

        # Construct message array: system prompt + user messages + end prompt
        all_messages = [
            {"role": "system", "content": CHAT_SYSTEM_PROMPT},
            *request.messages,
            {"role": "user", "content": CHAT_END_PROMPT}
        ]

        # Stream response
        try:
            stream = await self.client.chat.completions.create(
                model=model,
                messages=all_messages,
                stream=True
            )
            
            async for chunk in stream:
                if chunk.choices[0].delta.content:
                    yield chunk.choices[0].delta.content

        except Exception as e:
            logger.error(f"LLM streaming error: {e}")
            # Don't raise HTTPException here - streaming already started
            # Just log the error and end the stream
            yield f"\n[Error: Failed to generate complete response]"

