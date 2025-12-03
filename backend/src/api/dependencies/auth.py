"""
Authentication dependencies for FastAPI endpoints.
Enhanced with caching, admin impersonation, and care provider validation.
"""
import hashlib
import logging
from typing import Optional
from uuid import UUID

from authlib.jose import JoseError, jwt
from fastapi import Depends, Header, HTTPException, status
from fastapi.security import HTTPAuthorizationCredentials, HTTPBearer
from supabase import Client

from src.config.settings import get_settings

logger = logging.getLogger(__name__)

# HTTP Bearer token scheme
security = HTTPBearer()


# ============================================================================
# JWT Decoding & Validation
# ============================================================================


def decode_jwt_local(token: str) -> Optional[dict]:
    """
    Decode and validate Supabase JWT token locally (fast path).

    Returns claims if valid, None if invalid or expired.
    Avoids hitting Supabase API for every request.
    """
    settings = get_settings()

    if not token or not settings.supabase_jwt_secret:
        return None

    try:
        claims = jwt.decode(token, settings.supabase_jwt_secret)
        
        try:
            claims.validate(
                leeway=120,  # Allow 2 minutes clock skew
                claims_options={"aud": {"essential": False}},
            )
            return dict(claims)
        except Exception as e:
            logger.debug(f"JWT validation failed: {e}")
            return None
    except JoseError as e:
        logger.debug(f"JWT decode failed: {e}")
        return None
    except Exception as e:
        logger.warning(f"Unexpected error decoding JWT: {e}")
        return None
