"""
Supabase client management and dependency injection.
Provides schema-specific clients and dependency injection for FastAPI.
"""
from functools import lru_cache
from typing import Annotated

from fastapi import Depends
from supabase import Client, create_client
from supabase.client import ClientOptions

from src.config.settings import Settings, get_settings


@lru_cache(maxsize=None)
def get_supabase_for(
    postgrest_client_timeout: int = 60,
    storage_client_timeout: int = 60,
    schema: str = "public",
) -> Client:
    """
    Returns a Supabase client with the specified options.

    Args:
        postgrest_client_timeout (int): Timeout for PostgREST client in seconds.
        storage_client_timeout (int): Timeout for storage client in seconds.
        schema (str): The Postgres schema to use (defaults to "public").
                     Note: This is stored for reference but the Supabase Python client
                     requires explicit schema in table names for non-public schemas.

    Returns:
        Client: Configured Supabase client.

    Raises:
        ValueError: If SUPABASE_URL or SUPABASE_SERVICE_ROLE_KEY environment variables are not set.
    """
    settings = get_settings()

    if not settings.supabase_url or not settings.supabase_service_role_key:
        raise ValueError(
            "SUPABASE_URL and SUPABASE_SERVICE_ROLE_KEY must be set in environment variables"
        )

    # Note: The Supabase Python client's ClientOptions schema parameter doesn't actually
    # work as expected. We need to explicitly specify schema in table names.
    # See: https://github.com/supabase-community/supabase-py/issues
    options = ClientOptions(
        postgrest_client_timeout=postgrest_client_timeout,
        storage_client_timeout=storage_client_timeout,
        schema=schema,
    )
    client = create_client(
        settings.supabase_url, settings.supabase_service_role_key, options=options
    )
    # Store schema for later reference
    client._schema = schema  # Store schema on client for repository use
    return client


# Schema-specific client functions
def get_supabase() -> Client:
    """Get default Supabase client (public schema)."""
    return get_supabase_for(schema="public")


def get_supabase_public() -> Client:
    """Get Supabase client for public schema."""
    return get_supabase_for(schema="public")

def get_supabase_user_scoped(
    access_token: str,
    postgrest_client_timeout: int = 60,
    storage_client_timeout: int = 60,
) -> Client:
    """
    Get Supabase client with user-scoped access token for RLS.

    Args:
        access_token: User's JWT access token
        postgrest_client_timeout: Timeout for PostgREST client
        storage_client_timeout: Timeout for storage client

    Returns:
        Client: User-scoped Supabase client
    """
    settings = get_settings()

    client = create_client(
        settings.supabase_url,
        settings.supabase_key,  # Use anon key for user-scoped access
        options=ClientOptions(
            postgrest_client_timeout=postgrest_client_timeout,
            storage_client_timeout=storage_client_timeout,
            schema="public",
        ),
    )
    client.postgrest.auth(access_token)
    return client


# Type aliases for dependency injection
SupabaseClient = Annotated[Client, Depends(get_supabase)]
SupabasePublicClient = Annotated[Client, Depends(get_supabase_public)]
