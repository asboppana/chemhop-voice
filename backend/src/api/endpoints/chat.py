"""
LLM chat endpoints.

Provides health chat functionality with streaming responses.
"""
from fastapi import APIRouter, Depends, status
from fastapi.responses import StreamingResponse

from src.api.models import ErrorResponse, ChatRequest
from src.controllers.chat_controller import LLMController

# ============================================================================
# Dependency Injection
# ============================================================================


def get_llm_controller() -> LLMController:
    """Dependency injection for LLMController."""
    return LLMController()


# ============================================================================
# Router
# ============================================================================

router = APIRouter()


# ============================================================================
# Endpoints
# ============================================================================


@router.post(
    "/chat/general",
    status_code=status.HTTP_200_OK,
    responses={
        400: {"model": ErrorResponse, "description": "Invalid request"},
        401: {"model": ErrorResponse, "description": "Unauthorized"},
        403: {"model": ErrorResponse, "description": "Access denied"},
        500: {"model": ErrorResponse, "description": "Internal server error"},
    },
)
async def health_chat(
    request: ChatRequest,
    controller: LLMController = Depends(get_llm_controller),
) -> StreamingResponse:
    """
    Health chat endpoint with streaming LLM response.
    
    Takes health information as messages and returns an LLM-generated response
    using OpenAI API. Expects OPENAI_API_KEY environment variable.
    
    Requires authentication. Streams response progressively as text/plain.
    """
    # Validate request BEFORE creating StreamingResponse
    # This ensures validation errors (400) are raised before streaming starts
    controller._validate_request(request)
    
    return StreamingResponse(
        controller.health_chat_stream(request),
        media_type="text/plain",
    )

