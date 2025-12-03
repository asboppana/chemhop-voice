"""
Health check endpoints.
Simple endpoints for monitoring application health and status.
"""
from datetime import datetime

from fastapi import APIRouter, Depends
from pydantic import BaseModel

from src.config.settings import Settings, get_settings

router = APIRouter(tags=["health"])


class HealthResponse(BaseModel):
    """Health check response model."""

    status: str
    timestamp: datetime
    version: str
    environment: str


class DetailedHealthResponse(BaseModel):
    """Detailed health check response model."""

    status: str
    timestamp: datetime
    version: str
    environment: str
    services: dict
    database: dict


@router.get("/health", response_model=HealthResponse)
async def health_check(settings: Settings = Depends(get_settings)):
    """Basic health check endpoint."""
    return HealthResponse(
        status="ok",
        timestamp=datetime.utcnow(),
        version="2.0.0",
        environment=settings.environment,
    )


@router.get("/health/detailed", response_model=DetailedHealthResponse)
async def detailed_health_check(settings: Settings = Depends(get_settings)):
    """Detailed health check with service status."""
    # In production, these would be actual health checks
    services = {
        "email": "ok" if settings.resend_api_key else "not_configured",
        "redis": "ok",  # Would check Redis connection
        "twilio": "ok" if settings.twilio_account_sid else "not_configured",
        "stripe": "ok" if settings.stripe_secret_key else "not_configured",
    }

    database = {
        "status": "ok",  # Would check database connection
        "host": settings.supabase_host,
        "port": settings.supabase_db_port,
    }

    return DetailedHealthResponse(
        status="ok",
        timestamp=datetime.utcnow(),
        version="2.0.0",
        environment=settings.environment,
        services=services,
        database=database,
    )
