"""
Sedona Health Platform - Backend V2
Clean Architecture FastAPI Application with Supabase
"""
# Import realtime compatibility shim FIRST, before any supabase imports
# This patches realtime to add back AuthorizationError and NotConnectedError
# that were removed in realtime 2.24.0 but are still imported by supabase 2.24.0
import src.utils.realtime_compat  # noqa: F401

import os
import logging
from contextlib import asynccontextmanager
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles

from src.config.settings import get_settings
from src.api.routers import api_router
from src.middleware.request_logging import RequestLoggingMiddleware
from src.middleware.error_handling import ErrorHandlingMiddleware
from src.services.elevenlabs import create_elevenlabs_agent
from src.services.mcp.mcp_server import mcp


@asynccontextmanager
async def lifespan(app: FastAPI):
    """Application lifespan context manager."""
    # Startup
    logging.info("Starting Sedona Health API V2 with Supabase")
    
    # Validate Supabase configuration
    settings = get_settings()
    if not settings.supabase_url or not settings.supabase_service_role_key:
        logging.error("Supabase configuration missing! Check SUPABASE_URL and SUPABASE_SERVICE_ROLE_KEY")
    else:
        logging.info(f"Supabase configured for environment: {settings.environment}")
    
    # Initialize voice agent service
    try:
        app.state.voice_agent_id = create_elevenlabs_agent("drugdiscoveryagent")
    except Exception as e:
        logging.error(f"Failed to create voice agent: {e}")
    
    # Initialize background jobs
    yield
    
    # Shutdown
    logging.info("Shutting down...")


def create_app() -> FastAPI:
    """Create and configure FastAPI application."""
    app = FastAPI(
        title="US Hacks API V1",
        description="Clean Architecture Backend for US Hacks Platform with Supabase",
        version="1.0.0",
        lifespan=lifespan,
        redirect_slashes=False
    )
    
    # Configure CORS
    app.add_middleware(
        CORSMiddleware,
        allow_origins=["*"],
        allow_credentials=True,
        allow_methods=["*"],
        allow_headers=["*"],
    )
    
    # Add custom middleware
    app.add_middleware(ErrorHandlingMiddleware)
    # Enable request logging in both dev and production
    app.add_middleware(RequestLoggingMiddleware)
    
    # Mount static files (public directory for assets like VCF cards)
    public_dir = os.path.join(os.path.dirname(__file__), "public")
    if os.path.exists(public_dir):
        app.mount("/public", StaticFiles(directory=public_dir), name="public")
    
    # Include API router
    app.include_router(api_router, prefix="/api/v1")

    # Mount the MCP server
    app.mount("/mcp", mcp.sse_app())
    
    return app


app = create_app()

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(
        "main:app",
        host="0.0.0.0",
        port=8000,
        reload=True
    )