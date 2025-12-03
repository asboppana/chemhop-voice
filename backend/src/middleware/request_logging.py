"""
Request logging middleware.
Logs HTTP requests and responses for monitoring and debugging.
"""
import json
import logging
import time
from typing import Callable

from fastapi import Request, Response
from starlette.middleware.base import BaseHTTPMiddleware

logger = logging.getLogger(__name__)


class RequestLoggingMiddleware(BaseHTTPMiddleware):
    """Middleware for logging HTTP requests and responses."""

    def __init__(self, app, ignore_paths: tuple = ()):
        super().__init__(app)
        self.ignore_paths = ignore_paths or (
            "/health",
            "/docs",
            "/openapi.json",
            "/redoc",
            "/favicon.ico",
        )

    async def dispatch(self, request: Request, call_next: Callable) -> Response:
        # Skip logging for ignored paths
        if request.url.path in self.ignore_paths:
            return await call_next(request)

        start_time = time.time()

        # Log request
        request_log = {
            "type": "request",
            "method": request.method,
            "path": request.url.path,
            "query_params": dict(request.query_params),
            "client_ip": self._get_client_ip(request),
            "user_agent": request.headers.get("user-agent"),
            "timestamp": start_time,
        }

        # Get user ID if available (from auth header)
        user_id = await self._extract_user_id(request)
        if user_id:
            request_log["user_id"] = user_id

        # Log request body for non-GET requests (for debugging)
        if request.method in ("POST", "PUT", "PATCH", "DELETE"):
            body = await self._get_request_body(request)
            if body:
                request_log["body"] = body

        logger.info(json.dumps(request_log))

        # Process request
        response = await call_next(request)

        # Calculate processing time
        process_time = time.time() - start_time

        # Log response
        response_log = {
            "type": "response",
            "method": request.method,
            "path": request.url.path,
            "status_code": response.status_code,
            "process_time_ms": round(process_time * 1000, 2),
            "client_ip": self._get_client_ip(request),
            "timestamp": time.time(),
        }

        if user_id:
            response_log["user_id"] = user_id

        # Add response size if available
        if hasattr(response, "body"):
            response_log["response_size_bytes"] = len(response.body)

        # Log level based on status code
        if response.status_code >= 500:
            logger.error(json.dumps(response_log))
        elif response.status_code >= 400:
            logger.warning(json.dumps(response_log))
        else:
            logger.info(json.dumps(response_log))

        # Add processing time header
        response.headers["X-Process-Time"] = str(process_time)

        return response

    def _get_client_ip(self, request: Request) -> str:
        """Extract client IP address from request."""
        # Check for forwarded headers (when behind proxy)
        forwarded = request.headers.get("x-forwarded-for")
        if forwarded:
            return forwarded.split(",")[0].strip()

        real_ip = request.headers.get("x-real-ip")
        if real_ip:
            return real_ip

        # Fallback to direct client
        if request.client:
            return request.client.host

        return "unknown"

    async def _extract_user_id(self, request: Request) -> str:
        """Extract user ID from authorization header if present."""
        try:
            auth_header = request.headers.get("authorization")
            if not auth_header or not auth_header.startswith("Bearer "):
                return None

            token = auth_header.split(" ")[1]

            # Import here to avoid circular imports
            from src.api.dependencies.auth import decode_jwt_local

            claims = decode_jwt_local(token)
            if claims:
                return claims.get("sub")

        except Exception:
            # Don't fail request if user extraction fails
            pass

        return None

    async def _get_request_body(self, request: Request) -> dict:
        """
        Safely extract and parse request body.
        Returns None if body cannot be read or parsed.
        """
        try:
            # Check if body has already been consumed
            if hasattr(request.state, "body"):
                body_bytes = request.state.body
            else:
                # Read body and cache it for downstream handlers
                body_bytes = await request.body()
                request.state.body = body_bytes

            if not body_bytes:
                return None

            # Try to parse as JSON
            body_str = body_bytes.decode("utf-8")
            body_data = json.loads(body_str)

            # Sanitize sensitive fields
            sensitive_fields = {"password", "token", "secret", "api_key", "credit_card"}
            if isinstance(body_data, dict):
                return {
                    k: "***REDACTED***" if k.lower() in sensitive_fields else v
                    for k, v in body_data.items()
                }

            return body_data

        except Exception:
            # Don't fail request if body parsing fails
            return None
