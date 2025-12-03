"""
Error handling middleware.
Centralizes error handling and response formatting for Supabase-based backend.
"""
import json
import logging
import traceback
from typing import Callable

from fastapi import Request, Response, status
from fastapi.responses import JSONResponse
from postgrest import APIError as PostgrestError
from pydantic import ValidationError
from starlette.middleware.base import BaseHTTPMiddleware

logger = logging.getLogger(__name__)


class ErrorHandlingMiddleware(BaseHTTPMiddleware):
    """Middleware for centralized error handling and logging."""

    async def _get_request_body(self, request: Request) -> dict:
        """
        Safely extract request body for error logging.
        """
        try:
            if hasattr(request.state, "body"):
                body_bytes = request.state.body
            else:
                body_bytes = await request.body()
                request.state.body = body_bytes

            if not body_bytes:
                return None

            body_str = body_bytes.decode("utf-8")
            return json.loads(body_str)
        except Exception:
            return None

    async def dispatch(self, request: Request, call_next: Callable) -> Response:
        try:
            response = await call_next(request)
            return response

        except ValidationError as e:
            # Get request body for context
            body = await self._get_request_body(request)
            
            logger.warning(
                "Validation error",
                extra={
                    "path": request.url.path,
                    "method": request.method,
                    "errors": e.errors(),
                    "request_body": body,
                },
            )
            return JSONResponse(
                status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
                content={
                    "error": "Validation Error",
                    "message": "Invalid input data",
                    "details": e.errors(),
                },
            )

        except PostgrestError as e:
            # Get request body for context
            body = await self._get_request_body(request)
            
            logger.error(
                "Supabase API error",
                extra={
                    "path": request.url.path,
                    "method": request.method,
                    "error": str(e),
                    "code": getattr(e, "code", None),
                    "request_body": body,
                },
                exc_info=True,
            )
            return JSONResponse(
                status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
                content={
                    "error": "Database Error",
                    "message": "A database error occurred. Please try again later.",
                },
            )

        except Exception as e:
            # Get request body for context
            body = await self._get_request_body(request)
            
            # Get full traceback for dev mode
            tb_str = traceback.format_exc()
            
            # Get settings to check if we're in dev mode
            from src.config.settings import get_settings

            try:
                settings = get_settings()
                is_production = settings.is_production
            except Exception:
                is_production = True  # Default to production mode for safety

            logger.error(
                "Unhandled error",
                extra={
                    "path": request.url.path,
                    "method": request.method,
                    "error": str(e),
                    "error_type": type(e).__name__,
                    "request_body": body,
                    "traceback": tb_str if not is_production else None,
                },
                exc_info=True,
            )

            # Don't expose internal errors in production
            if is_production:
                message = "An internal error occurred. Please try again later."
            else:
                # In dev mode, include full error details
                message = f"{type(e).__name__}: {str(e)}"

            response_content = {
                "error": "Internal Server Error",
                "message": message,
            }

            # Include traceback in dev mode
            if not is_production:
                response_content["traceback"] = tb_str

            return JSONResponse(
                status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
                content=response_content,
            )
