"""
Semantic Scholar API client for research endpoints.
"""
import os
from typing import Dict, Any, Optional, List
from aiohttp import ClientTimeout
import aiohttp

from src.utils.rate_limiter import RateLimiter

SEMANTIC_SCHOLAR_BASE_URL = "https://api.semanticscholar.org"
SEMANTIC_SCHOLAR_API_KEY = os.environ.get("SEMANTIC_SCHOLAR_API_KEY", "")
SEMANTIC_SCHOLAR_FIELDS = "paperId,corpusId,externalIds,url,title,abstract,venue,publicationVenue,year,referenceCount,citationCount,influentialCitationCount,isOpenAccess,openAccessPdf,fieldsOfStudy,s2FieldsOfStudy,publicationTypes,publicationDate,journal,citationStyles,authors"

# Rate limiter for Semantic Scholar API (1 request per second)
rate_limiter = RateLimiter(calls_per_second=1.0)


async def paper_search(
    session: aiohttp.ClientSession,
    query: str,
    limit: int = 25,
    year: Optional[str] = None,
    venue: Optional[str] = None,
    open_access: Optional[bool] = None,
    fields_of_study: Optional[List[str]] = None,
    publication_types: Optional[List[str]] = None,
    min_citation_count: Optional[int] = None,
) -> Dict[str, Any]:
    """
    Search for papers using Semantic Scholar API.
    
    Args:
        session: aiohttp session
        query: Search query string
        limit: Maximum number of results
        year: Publication year filter
        venue: Venue filter
        open_access: Filter for open access papers
        fields_of_study: List of fields of study
        publication_types: List of publication types
        min_citation_count: Minimum citation count
        
    Returns:
        Dictionary with 'data' key containing list of papers
    """
    url = f"{SEMANTIC_SCHOLAR_BASE_URL}/graph/v1/paper/search"
    
    params = {
        "query": query,
        "offset": 0,
        "limit": limit,
        "fields": SEMANTIC_SCHOLAR_FIELDS,
    }
    
    if year:
        params["year"] = year
    if venue:
        params["venue"] = venue
    if min_citation_count:
        params["minCitationCount"] = str(min_citation_count)
    if open_access is True:
        params["openAccessPdf"] = "openAccessPdf"
    if fields_of_study:
        params["fieldsOfStudy"] = ",".join(fields_of_study)
    if publication_types:
        params["publicationTypes"] = ",".join(publication_types)
    
    headers = {}
    if SEMANTIC_SCHOLAR_API_KEY:
        headers["x-api-key"] = SEMANTIC_SCHOLAR_API_KEY
    
    client_timeout = ClientTimeout(total=60)
    await rate_limiter.wait()
    
    async with session.get(url, params=params, headers=headers, timeout=client_timeout) as response:
        response.raise_for_status()
        return await response.json()

