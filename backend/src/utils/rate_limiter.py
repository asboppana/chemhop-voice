"""
Rate limiter utility for API rate limiting.
"""
import asyncio
import time
from typing import Optional


class RateLimiter:
    """
    Rate limiter for API calls.
    
    Ensures calls are spaced out to respect rate limits by waiting
    between requests if necessary.
    
    Example:
        limiter = RateLimiter(calls_per_second=1.0)
        await limiter.wait()  # Wait if needed before making request
        # Make API call here
    """
    
    def __init__(self, calls_per_second: float = 1.0):
        """
        Initialize rate limiter.
        
        Args:
            calls_per_second: Maximum number of calls allowed per second
        """
        self.calls_per_second = calls_per_second
        self.last_call: Optional[float] = None

    async def wait(self):
        """
        Wait if necessary to respect rate limit.
        
        Calculates time since last call and sleeps if needed to ensure
        we don't exceed the rate limit.
        """
        now = time.time()
        if self.last_call is not None:
            time_since_last_call = now - self.last_call
            min_interval = 1.0 / self.calls_per_second
            if time_since_last_call < min_interval:
                await asyncio.sleep(min_interval - time_since_last_call)
        self.last_call = time.time()

    def reset(self):
        """Reset the rate limiter (clear last call time)."""
        self.last_call = None

