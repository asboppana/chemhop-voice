from fastapi import APIRouter

from .endpoints import health
from .endpoints import chat

api_router = APIRouter()

# Include endpoint routers
# Health (no prefix)
api_router.include_router(health.router, prefix="", tags=["health"])
api_router.include_router(chat.router, prefix="", tags=["chat"])
