"""
Compatibility shim for realtime AuthorizationError and NotConnectedError.

Supabase 2.24.0 imports AuthorizationError and NotConnectedError from realtime,
but realtime 2.24.0 removed these exceptions. This module patches realtime
to add them back for compatibility.
"""
import realtime


# Define the missing exceptions if they don't exist
if not hasattr(realtime, 'AuthorizationError'):
    class AuthorizationError(Exception):
        """Raised when authorization fails."""
        pass
    
    realtime.AuthorizationError = AuthorizationError


if not hasattr(realtime, 'NotConnectedError'):
    class NotConnectedError(Exception):
        """Raised when not connected to realtime."""
        pass
    
    realtime.NotConnectedError = NotConnectedError

