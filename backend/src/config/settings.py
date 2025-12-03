"""
Application configuration and settings.
Centralized configuration management using Pydantic Settings.
"""
from functools import lru_cache
from typing import List, Optional

from pydantic import Field
from pydantic_settings import BaseSettings, SettingsConfigDict


class Settings(BaseSettings):
    """Application settings loaded from environment variables."""
    
    # Allow extra fields from .env files that aren't defined here
    # This is useful for variables used by other parts of the app (e.g., Supabase CLI, Logfire, etc.)
    model_config = SettingsConfigDict(
        extra="ignore",  # Ignore extra fields from .env files
        case_sensitive=False,
        env_file=[".env", ".env.local"],  # Load .env first, then .env.local (so .env.local overrides)
    )

    # Application settings
    app_name: str = "Sedona Health API V2"
    environment: str = Field(default="local", env="SYSTEM_ENVIRONMENT")
    debug: bool = Field(default=True, env="DEBUG")

    # Database settings
    supabase_url: str = Field(default="http://localhost:54321", env="SUPABASE_URL")
    supabase_key: str = Field(default="test_key", env="SUPABASE_KEY")
    supabase_service_role_key: str = Field(
        default="test_service_role_key", env="SUPABASE_SERVICE_ROLE_KEY"
    )
    supabase_host: str = Field(default="localhost", env="SUPABASE_HOST")
    supabase_db_port: str = Field(default="5432", env="SUPABASE_DB_PORT")

    # Authentication settings
    supabase_jwt_secret: str = Field(
        default="test_jwt_secret", env="SUPABASE_JWT_SECRET"
    )
    access_token_expire_minutes: int = Field(
        default=30, env="ACCESS_TOKEN_EXPIRE_MINUTES"
    )

    # CORS settings
    allowed_origins: Optional[List[str]] = None

    # Logging settings
    log_level: str = Field(default="INFO", env="LOG_LEVEL")
    enable_request_logging: bool = Field(default=False, env="ENABLE_REQUEST_LOGGING")

    # OpenAI settings
    openai_api_key: str = Field(default="", env="OPENAI_API_KEY")
    elevenlabs_api_key: str = Field(default="", env="ELEVENLABS_API_KEY")

    @property
    def is_production(self) -> bool:
        """Check if running in production environment."""
        return self.environment == "production"

    @property
    def is_local(self) -> bool:
        """Check if running in local development environment."""
        return self.environment == "local"


@lru_cache()
def get_settings() -> Settings:
    """Get cached application settings."""
    return Settings()
