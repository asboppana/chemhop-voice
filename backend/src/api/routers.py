from fastapi import APIRouter

from .endpoints import health
from .endpoints import molecule

api_router = APIRouter()

# Include endpoint routers
# Health (no prefix)
api_router.include_router(health.router, prefix="", tags=["health"])
api_router.include_router(molecule.router, prefix="", tags=["molecule"])
