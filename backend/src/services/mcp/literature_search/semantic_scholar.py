import os
import json
import argparse
import asyncio
from aiohttp import ClientTimeout
import aiohttp
import urllib3
from urllib.parse import urlencode
import time
from typing import List, Dict, Any, Union, Optional, AsyncGenerator
from enum import Enum

# Disable SSL warnings
urllib3.disable_warnings()

# Constants
SEMANTIC_SCHOLAR_BASE_URL = "https://api.semanticscholar.org"
SEMANTIC_SCHOLAR_API_KEY = os.environ.get("SEMANTIC_SCHOLAR_API_KEY", "")
if SEMANTIC_SCHOLAR_API_KEY:
    masked = SEMANTIC_SCHOLAR_API_KEY[:4] + "***" + SEMANTIC_SCHOLAR_API_KEY[-4:]
    print(f"✓ Semantic Scholar API key loaded: {masked}")
else:
    print("⚠ WARNING: SEMANTIC_SCHOLAR_API_KEY not found in environment. API requests will fail with 403 Forbidden.")
SEMANTIC_SCHOLAR_FIELDS = "paperId,corpusId,externalIds,url,title,abstract,venue,publicationVenue,year,referenceCount,citationCount,influentialCitationCount,isOpenAccess,openAccessPdf,fieldsOfStudy,s2FieldsOfStudy,publicationTypes,publicationDate,journal,citationStyles,authors,citations,references,embedding,tldr"
SEMANTIC_SCHOLAR_RECOMMENDATIONS_FIELDS = "paperId,corpusId,externalIds,url,title,abstract,venue,publicationVenue,year,referenceCount,citationCount,influentialCitationCount,isOpenAccess,openAccessPdf,fieldsOfStudy,s2FieldsOfStudy,publicationTypes,publicationDate,journal,citationStyles,authors"
SEMANTIC_SCHOLAR_AUTHOR_FIELDS = "authorId,externalIds,url,name,affiliations,homepage,paperCount,citationCount,hIndex"
SEMANTIC_SCHOLAR_CONTEXT_FIELDS = "contexts,intents,contextsWithIntent,isInfluential"


class PublicationTypes(str, Enum):
    Review: str
    JournalArticle: str
    CaseReport: str
    ClinicalTrial: str
    Conference: str
    Dataset: str
    Editorial: str
    LettersAndComments: str
    MetaAnalysis: str
    News: str
    Study: str
    Book: str
    BookSection: str


class FieldOfStudy(str, Enum):
    ComputerScience: str = "Computer Science"
    Medicine: str = "Medicine"
    Chemistry: str = "Chemistry"
    Biology: str = "Biology"
    MaterialsScience: str = "Materials Science"
    Physics: str = "Physics"
    Geology: str = "Geology"
    Psychology: str = "Psychology"
    Art: str = "Art"
    History: str = "History"
    Geography: str = "Geography"
    Sociology: str = "Sociology"
    Business: str = "Business"
    PoliticalScience: str = "Political Science"
    Economics: str = "Economics"
    Philosophy: str = "Philosophy"
    Mathematics: str = "Mathematics"
    Engineering: str = "Engineering"
    EnvironmentalScience: str = "Environmental Science"
    AgriculturalAndFoodSciences: str = "Agricultural and Food Sciences"
    Education: str = "Education"
    Law: str = "Law"
    Linguistics: str = "Linguistics"


class Paper(str, Enum):
    paperId = "paperId"
    corpusId = "corpusId"
    externalIds = "externalIds"
    url = "url"
    title = "title"
    abstract = "abstract"
    venue = "venue"
    publicationVenue = "publicationVenue"
    year = "year"
    referenceCount = "referenceCount"
    citationCount = "citationCount"
    isInfluential = "isInfluential"
    influentialCitationCount = "influentialCitationCount"
    isOpenAccess = "isOpenAccess"
    openAccessPdf = "openAccessPdf"
    fieldsOfStudy = "fieldsOfStudy"
    s2FieldsOfStudy = "s2FieldsOfStudy"
    publicationTypes = "publicationTypes"
    publicationDate = "publicationDate"
    journal = "journal"
    citationStyles = "citationStyles"
    authors = "authors"
    citations = "citations"
    references = "references"
    embedding = "embedding"
    tldr = "tldr"


class PaperID(str, Enum):
    """
    <sha> - a Semantic Scholar ID, e.g. 649def34f8be52c8b66281af98ae884c09aef38b
    CorpusId:<id> - a Semantic Scholar numerical ID, e.g. CorpusId:215416146
    DOI:<doi> - a Digital Object Identifier, e.g. DOI:10.18653/v1/N18-3011
    ARXIV:<id> - arXiv.rg, e.g. ARXIV:2106.15928
    MAG:<id> - Microsoft Academic Graph, e.g. MAG:112218234
    ACL:<id> - Association for Computational Linguistics, e.g. ACL:W12-3903
    PMID:<id> - PubMed/Medline, e.g. PMID:19872477
    PMCID:<id> - PubMed Central, e.g. PMCID:2323736
    URL:<url> - URL from one of the sites listed below
        semanticscholar.org
        arxiv.org
        aclweb.org
        acm.org
        biorxiv.org
    """

    SHA = "sha"
    CORPUS_ID = "CorpusId"
    DOI = "DOI"
    ARXIV = "ARXIV"
    MAG = "MAG"
    ACL = "ACL"
    PMID = "PMID"
    PMCID = "PMCID"
    URL = "URL"


class Author(str, Enum):
    authorId = "authorId"
    externalIds = "externalIds"
    url = "url"
    name = "name"
    affiliations = "affiliations"
    homepage = "homepage"
    paperCount = "paperCount"
    citationCount = "citationCount"
    hIndex = "hIndex"


class Embedding(str, Enum):
    model = "model"
    vector = "vector"


class Tldr(str, Enum):
    model = "model"
    text = "text"


class Context(str, Enum):
    contexts = "contexts"
    intents = "intents"
    contexts_with_intent = "contextsWithIntent"
    is_influential = "isInfluential"
    citing_paper = "citingPaper"
    cited_paper = "citedPaper"


class RateLimiter:
    """
    Rate limiting for Semantic Scholar API.
    
    CRITICAL: 1 request per second MAX for:
    - /paper/batch
    - /paper/search
    - /recommendations
    
    10 requests/second for all other calls.
    
    We use 0.9 calls/sec (~1.11s between calls) as a safety buffer to ensure
    we NEVER exceed the 1/sec limit, which can result in 403 Forbidden errors.
    """

    def __init__(self, calls_per_second: float = 1.0):
        self.calls_per_second = calls_per_second
        self.last_call = 0

    async def wait(self):
        now = time.time()
        time_since_last_call = now - self.last_call
        if time_since_last_call < 1 / self.calls_per_second:
            await asyncio.sleep(1 / self.calls_per_second - time_since_last_call)
        self.last_call = time.time()


rate_limiter = RateLimiter(calls_per_second=0.9)  # 0.9 calls/sec = ~1.11s between calls (safety buffer)


# Helper functions
async def make_request(
    session: aiohttp.ClientSession,
    url: str,
    params: Dict[str, Any] = None,
    method: str = "GET",
    data: Dict[str, Any] = None,
    timeout: int = 600,
    retry_count: int = 0,
) -> Dict[str, Any]:
    headers = {"x-api-key": SEMANTIC_SCHOLAR_API_KEY} if SEMANTIC_SCHOLAR_API_KEY else {}
    client_timeout = ClientTimeout(total=timeout)
    await rate_limiter.wait()  # Wait before making the request
    print(f"Making {method} request to {url} with params: {params}")

    try:
        if method == "GET":
            async with session.get(url, params=params, headers=headers, timeout=client_timeout) as response:
                print(response.status)
                if response.status >= 400:
                    try:
                        body = await response.text()
                        print(f"Response body: {body[:500]}")
                    except Exception:
                        pass
                    response.raise_for_status()
                return await response.json()
        elif method == "POST":
            headers["Content-Type"] = "application/json"
            async with session.post(url, params=params, headers=headers, json=data, timeout=client_timeout) as response:
                if response.status >= 400:
                    try:
                        body = await response.text()
                        print(f"Response body: {body[:500]}")
                    except Exception:
                        pass
                    response.raise_for_status()
                return await response.json()
    except asyncio.TimeoutError:
        print(f"Request timed out after {timeout} seconds")
        raise
    except aiohttp.ClientResponseError as e:
        print(f"HTTP error occurred: {e.status} {e.message}")
        if e.status == 429 and retry_count < 3:  # Limit retries to prevent infinite recursion
            wait_time = 5 * (retry_count + 1)  # Exponential backoff
            print(f"Rate limit exceeded. Waiting for {wait_time} seconds before retrying...")
            await asyncio.sleep(wait_time)
            return await make_request(session, url, params, method, data, timeout, retry_count + 1)
        raise
    except aiohttp.ClientError as e:
        print(f"An error occurred: {str(e)}")
        raise


async def paginate(session: aiohttp.ClientSession, url: str, params: Dict[str, Any], data_key: str = "data", timeout: int = 300) -> AsyncGenerator[Dict[str, Any], None]:
    while True:
        response = await make_request(session, url, params, timeout=timeout)
        for item in response.get(data_key, []):
            yield item

        next_offset = response.get("next")
        if next_offset is None:
            break

        params["offset"] = next_offset


async def _prepare_search_params(
    query: str = None,
    paper_id: str = None,
    offset: int = 0,
    limit: int = 100,
    fields: Optional[str] = None,
    year: Optional[Union[int, str]] = None,
    venue: Optional[str] = None,
    open_access: Optional[bool] = None,
    fields_of_study: Optional[List[str]] = None,
    min_citation_count: Optional[int] = None,
    publication_types: Optional[List[str]] = None,
    type: str = "paper",
) -> Dict[str, Any]:
    params = {
        "offset": offset,
        "limit": limit,
        "fields": fields or (SEMANTIC_SCHOLAR_FIELDS if type == "paper" else SEMANTIC_SCHOLAR_RECOMMENDATIONS_FIELDS),
    }
    if query:
        params["query"] = query
    if paper_id:
        params["paperId"] = paper_id
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
    return params


"""
Paper Search Functions
"""


async def paper_autocomplete(session: aiohttp.ClientSession, query: str) -> List[Dict[str, Any]]:
    url = f"{SEMANTIC_SCHOLAR_BASE_URL}/graph/v1/paper/autocomplete"
    params = {"query": query}
    response = await make_request(session, url, params)
    return response.get("matches", [])


async def paper_batch(session: aiohttp.ClientSession, paper_ids: List[str], fields: str = None) -> List[Dict[str, Any]]:
    url = f"{SEMANTIC_SCHOLAR_BASE_URL}/graph/v1/paper/batch"  # 500 paper IDs at a time
    params = {"fields": fields or SEMANTIC_SCHOLAR_FIELDS}
    data = {"ids": paper_ids}
    return await make_request(session, url, params, method="POST", data=data)


async def paper_search(session: aiohttp.ClientSession, **kwargs) -> Dict[str, Any]:
    url = f"{SEMANTIC_SCHOLAR_BASE_URL}/graph/v1/paper/search"
    params = await _prepare_search_params(**kwargs)
    return await make_request(session, url, params)


async def search_papers_generator(session: aiohttp.ClientSession, **kwargs) -> AsyncGenerator[Dict[str, Any], None]:
    url = f"{SEMANTIC_SCHOLAR_BASE_URL}/graph/v1/paper/search"
    params = await _prepare_search_params(**kwargs)
    async for paper in paginate(session, url, params):
        yield paper


async def paper_search_bulk(session: aiohttp.ClientSession, **kwargs) -> Dict[str, Any]:
    url = f"{SEMANTIC_SCHOLAR_BASE_URL}/graph/v1/paper/search/bulk"
    params = await _prepare_search_params(**kwargs)
    return await make_request(session, url, params)


async def paper_search_match(session: aiohttp.ClientSession, **kwargs) -> Dict[str, Any]:
    url = f"{SEMANTIC_SCHOLAR_BASE_URL}/graph/v1/paper/search/match"
    params = await _prepare_search_params(**kwargs)
    return await make_request(session, url, params)


"""
Paper Details
"""


async def get_paper_details(session: aiohttp.ClientSession, paper_id: str, fields: str = None) -> Dict[str, Any]:
    url = f"{SEMANTIC_SCHOLAR_BASE_URL}/graph/v1/paper/{paper_id}"
    params = {"fields": fields or SEMANTIC_SCHOLAR_FIELDS}
    return await make_request(session, url, params)


async def get_paper_authors(session: aiohttp.ClientSession, paper_id: str, offset: int = 0, limit: int = 100, fields: str = None) -> Dict[str, Any]:
    url = f"{SEMANTIC_SCHOLAR_BASE_URL}/graph/v1/paper/{paper_id}/authors"
    params = {"offset": offset, "limit": limit, "fields": fields or SEMANTIC_SCHOLAR_AUTHOR_FIELDS}
    return await make_request(session, url, params)


async def get_paper_citations(session: aiohttp.ClientSession, paper_id: str, offset: int = 0, limit: int = 100, fields: str = None) -> Dict[str, Any]:
    url = f"{SEMANTIC_SCHOLAR_BASE_URL}/graph/v1/paper/{paper_id}/citations"  # papers citing this paper
    fields = fields or SEMANTIC_SCHOLAR_CONTEXT_FIELDS + ",citingPaper"  # Include details about the citing paper
    params = {"offset": offset, "limit": limit, "fields": fields}
    return await make_request(session, url, params)


async def get_paper_references(session: aiohttp.ClientSession, paper_id: str, offset: int = 0, limit: int = 100, fields: str = None) -> Dict[str, Any]:
    url = f"{SEMANTIC_SCHOLAR_BASE_URL}/graph/v1/paper/{paper_id}/references"  # papers cited by this paper
    fields = fields or SEMANTIC_SCHOLAR_CONTEXT_FIELDS + ",citedPaper"  # Include details about the cited paper
    params = {"offset": offset, "limit": limit, "fields": fields}
    return await make_request(session, url, params)


"""
Author Functions
"""


async def author_batch(session: aiohttp.ClientSession, author_ids: List[str], fields: str = None) -> List[Dict[str, Any]]:
    url = f"{SEMANTIC_SCHOLAR_BASE_URL}/graph/v1/author/batch"
    params = {"fields": fields or SEMANTIC_SCHOLAR_AUTHOR_FIELDS}
    data = {"ids": author_ids}
    return await make_request(session, url, params, method="POST", data=data)


async def author_search(session: aiohttp.ClientSession, query: str, offset: int = 0, limit: int = 100, fields: str = None) -> Dict[str, Any]:
    url = f"{SEMANTIC_SCHOLAR_BASE_URL}/graph/v1/author/search"
    params = {"query": query, "offset": offset, "limit": limit, "fields": fields or SEMANTIC_SCHOLAR_AUTHOR_FIELDS}
    return await make_request(session, url, params)


async def get_author_details(session: aiohttp.ClientSession, author_id: str, fields: str = None) -> Dict[str, Any]:
    url = f"{SEMANTIC_SCHOLAR_BASE_URL}/graph/v1/author/{author_id}"
    params = {"fields": fields or SEMANTIC_SCHOLAR_AUTHOR_FIELDS}
    return await make_request(session, url, params)


async def get_author_papers(session: aiohttp.ClientSession, author_id: str, offset: int = 0, limit: int = 100, fields: str = None) -> Dict[str, Any]:
    url = f"{SEMANTIC_SCHOLAR_BASE_URL}/graph/v1/author/{author_id}/papers"
    params = {"offset": offset, "limit": limit, "fields": fields or SEMANTIC_SCHOLAR_FIELDS}
    return await make_request(session, url, params)


"""
Recommendation Functions
"""


async def get_paper_recommendations(session: aiohttp.ClientSession, positive_paper_ids: List[str], negative_paper_ids: List[str] = None, **kwargs) -> List[Dict[str, Any]]:
    url = f"{SEMANTIC_SCHOLAR_BASE_URL}/recommendations/v1/papers"
    params = await _prepare_search_params(**kwargs, type="recommendations")
    data = {"positivePaperIds": positive_paper_ids, "negativePaperIds": negative_paper_ids or []}
    return await make_request(session, url, params, method="POST", data=data)


async def get_single_paper_recommendations(session: aiohttp.ClientSession, paper_id: str, from_pool: str = "recent", **kwargs) -> List[Dict[str, Any]]:
    url = f"{SEMANTIC_SCHOLAR_BASE_URL}/recommendations/v1/papers/forpaper/{paper_id}"
    params = await _prepare_search_params(**kwargs, type="recommendations")
    params["from"] = from_pool  # all-cs
    return await make_request(session, url, params)


"""
Dataset Functions
"""


async def get_available_releases(session: aiohttp.ClientSession) -> List[str]:
    url = f"{SEMANTIC_SCHOLAR_BASE_URL}/datasets/v1/release/"
    return await make_request(session, url)


async def get_release_info(session: aiohttp.ClientSession, release_id: str) -> Dict[str, Any]:
    url = f"{SEMANTIC_SCHOLAR_BASE_URL}/datasets/v1/release/{release_id}"
    return await make_request(session, url)


async def get_dataset_download_links(session: aiohttp.ClientSession, release_id: str, dataset_name: str) -> Dict[str, Any]:
    url = f"{SEMANTIC_SCHOLAR_BASE_URL}/datasets/v1/release/{release_id}/dataset/{dataset_name}"
    return await make_request(session, url)


async def get_dataset_diffs(session: aiohttp.ClientSession, dataset_name: str, start_release_id: str, end_release_id: str) -> Dict[str, Any]:
    url = f"{SEMANTIC_SCHOLAR_BASE_URL}/datasets/v1/diffs/{start_release_id}/to/{end_release_id}/{dataset_name}"
    return await make_request(session, url)


async def update_dataset(session: aiohttp.ClientSession, dataset_name: str, current_release: str, target_release: str):
    diffs = await get_dataset_diffs(session, dataset_name, current_release, target_release)
    print(f"Found {len(diffs)} differences between releases {current_release} and {target_release}")


"""
Main Function
"""


async def main(args):
    BLUE = "\033[94m"
    RED = "\033[91m"
    RESET = "\033[0m"
    print(f"{BLUE}Args: {args}{RESET}")
    search_params = {
        "query": args.query,
        "fields": args.fields,
        "year": args.year,
        "venue": args.venue,
        "min_citation_count": args.min_citation_count,
        "open_access": args.open_access,
        "fields_of_study": args.fields_of_study,
        "publication_types": args.publication_types,
    }

    async with aiohttp.ClientSession() as session:
        if args.action == "autocomplete":
            completions = await paper_autocomplete(session, args.query)
            print(json.dumps(completions, indent=2))

        elif args.action == "search":
            if args.paginate:
                async for paper in search_papers_generator(session, **search_params):
                    print(json.dumps(paper, indent=2))
            else:
                print(f"{RED}Search params: {search_params}{RESET}")
                papers = await paper_search(session, limit=args.limit, **search_params)
                print(json.dumps(papers, indent=2))
        elif args.action == "search_bulk":
            papers = await paper_search_bulk(session, limit=args.limit, **search_params)
            print(json.dumps(papers, indent=2))
        elif args.action == "match":  # deprecated in favor of recommendations
            papers = await paper_search_match(session, limit=args.limit, **search_params)
            print(json.dumps(papers, indent=2))
        elif args.action == "batch":
            papers = await paper_batch(session, args.paper_ids, args.fields)
            print(json.dumps(papers, indent=2))

        elif args.action == "details":
            paper = await get_paper_details(session, args.paper_id, args.fields)
            print(json.dumps(paper, indent=2))
        elif args.action == "citing_details":
            citations = await get_paper_citations(session, args.paper_id, limit=args.limit, fields=args.fields)
            print(json.dumps(citations, indent=2))
        elif args.action == "cited_details":
            references = await get_paper_references(session, args.paper_id, limit=args.limit, fields=args.fields)
            print(json.dumps(references, indent=2))

        elif args.action == "author_search":
            authors = await author_search(session, args.query, fields=args.fields, limit=args.limit)
            print(json.dumps(authors, indent=2))
        elif args.action == "author_details":
            author = await get_author_details(session, args.author_id, args.fields)
            print(json.dumps(author, indent=2))

        elif args.action == "recommend":
            if args.paper_id:
                recommendations = await get_single_paper_recommendations(session, args.paper_id, args.from_pool, **search_params)
            else:
                recommendations = await get_paper_recommendations(session, args.pos_paper_ids, args.neg_paper_ids, **search_params)
            print(json.dumps(recommendations, indent=2))

        elif args.action == "list_releases":
            releases = await get_available_releases(session)
            print(json.dumps(releases, indent=2))
        elif args.action == "release_info":
            release_info = await get_release_info(session, args.release_id)
            print(json.dumps(release_info, indent=2))
        elif args.action == "update_dataset":
            await update_dataset(session, args.dataset_name, args.current_release, args.target_release)


if __name__ == "__main__":
    """
    python3 semantic_scholar.py search --query 'high intensity interval training HIIT LDL cholesterol reduction' --fields-of-study Medicine --publication-types MetaAnalysis ClinicalTrial --min-citation-count 10 --open-access --limit 50
    """
    parser = argparse.ArgumentParser(description="Semantic Scholar API Client")
    parser.add_argument(
        "action",
        choices=[
            "autocomplete",
            "search",
            "search_bulk",
            "match",
            "batch",
            "details",
            "citing_details",
            "cited_details",
            "author_search",
            "author_details",
            "recommend",
            "list_releases",
            "release_info",
            "update_dataset",
        ],
        help="Action to perform",
    )

    # For all: what to return from a paper or author object
    parser.add_argument("--fields", help="Comma-separated list of fields to return")
    # How many results to return
    parser.add_argument("--limit", type=int, default=100, help="Limit the number of results")
    parser.add_argument("--paginate", action="store_true", help="Use pagination for search results")

    # For paper search
    parser.add_argument("--query", help="Search query for papers")
    parser.add_argument("--year", help="Publication year or range (e.g., 2020 or 2018-2021 or 2019-03-05:)")
    parser.add_argument("--venue", help="Publication venue")
    parser.add_argument("--min-citation-count", type=int, help="Minimum citation count")
    parser.add_argument("--open-access", action="store_true", help="Filter for open access papers")
    parser.add_argument("--fields-of-study", nargs="+", help="Fields of study to filter by")
    parser.add_argument("--publication-types", nargs="+", help="Publication types to filter by")

    # For bulk
    parser.add_argument("--paper-ids", nargs="+", help="Paper IDs for batch lookup")

    # For paper and author details
    parser.add_argument("--paper-id", help="Paper ID")
    parser.add_argument("--author-id", help="Author ID")

    # For recommendations
    parser.add_argument("--pos-paper-ids", nargs="+", help="Positive paper IDs for recommendations")
    parser.add_argument("--neg-paper-ids", nargs="+", help="Negative paper IDs for recommendations")
    parser.add_argument("--from-pool", choices=["recent", "all-cs"], default="recent", help="Which pool of papers to recommend from")

    # For bulk download
    parser.add_argument("--release-id", help="Release ID for dataset operations")
    parser.add_argument("--dataset-name", help="Name of the dataset")
    parser.add_argument("--current-release", help="Current release ID for dataset update")
    parser.add_argument("--target-release", help="Target release ID for dataset update")

    args = parser.parse_args()
    asyncio.run(main(args))
