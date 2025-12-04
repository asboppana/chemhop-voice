from fastapi import APIRouter

from .endpoints import health
from .endpoints import molecule
from .endpoints import voice_agent
from .endpoints import admet

api_router = APIRouter()

# Include endpoint routers
# Health (no prefix)
api_router.include_router(health.router, prefix="", tags=["health"])
api_router.include_router(molecule.router, prefix="", tags=["molecule"])
api_router.include_router(voice_agent.router, prefix="", tags=["voice_agent"])
api_router.include_router(admet.router, prefix="", tags=["admet"])
