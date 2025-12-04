"""
MCP Server for Drug Discovery Tools.

This server exposes chemistry and drug discovery tools via the Model Context Protocol (MCP).
It provides access to:
- Smart Chemist: Annotate molecular structures with functional groups
- Ring Scanner: Find bio-isosteric ring replacements

Note: ADMET AI predictions are now available via a separate MCP server (admet-service).

Run with: python -m src.services.mcp.mcp_server
"""

import sys
import logging
from pathlib import Path
from typing import Any, Dict

from rdkit import Chem

from typing import Optional

import requests
from fastmcp import FastMCP

# Configure logging with colors
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%H:%M:%S'
)
logger = logging.getLogger(__name__)

# Color codes for terminal output
_COLOR_GREEN = "\033[92m"
_COLOR_BLUE = "\033[94m"
_COLOR_YELLOW = "\033[93m"
_COLOR_RESET = "\033[0m"

# Add the backend directory to sys.path for imports
backend_dir = Path(__file__).resolve().parent.parent.parent.parent
if str(backend_dir) not in sys.path:
    sys.path.insert(0, str(backend_dir))

# Import our tools
from src.services.mcp.smart_chemist.tools.smart_chemist import (  # noqa: E402
    SmartChemist,
    convert_string_input_to_smiles,
)
from src.services.mcp.patent_search.surechembl_tool import SureChEMBLPatentTool
from src.services.mcp.ring_scan.tools.bioisostere_scanner import get_scanner  # noqa: E402

# Initialize FastMCP server
mcp = FastMCP("drugdiscovery_mcp")

# Initialize tool instances
smart_chemist = SmartChemist()
patent_tool = SureChEMBLPatentTool()
ring_scanner = get_scanner()


@mcp.tool()
def annotate_molecule(smiles: str) -> Dict[str, Any]:
    """
    Annotate a molecule with functional groups and structural patterns.
    
    This tool analyzes a molecule's structure and identifies known functional groups,
    ring systems, and chemical patterns using a comprehensive database of SMARTS patterns.
    It's useful for:
    - Understanding molecular structure and composition
    - Identifying key functional groups for drug-likeness
    - Finding specific chemical moieties
    - Educational purposes to learn about molecular features
    
    Args:
        smiles: SMILES string representation of the molecule
                Examples: "CCO" (ethanol), "c1ccccc1" (benzene), 
                         "CC(=O)Oc1ccccc1C(=O)O" (aspirin)
    
    Returns:
        Dictionary containing:
        - name: Molecule name (if available in the input)
        - svg: SVG image representation of the molecule
        - smiles: The canonical SMILES string
        - matches: List of matched patterns, each containing:
            - atom_indices: Indices of atoms matching the pattern
            - trivial_name: Pattern information including:
                - name: Common name of the functional group
                - smarts: SMARTS pattern
                - group: Category (e.g., "cyclic", "acyclic", "functional_group")
                - bonds: Number of bonds in the pattern
                - hierarchy: Hierarchical relationship to other patterns
    
    Raises:
        ValueError: If the SMILES string is invalid
    
    Example:
        >>> annotate_molecule("CCO")
        {
            "name": "ethanol",
            "svg": "<svg>...</svg>",
            "smiles": "CCO",
            "matches": [
                {
                    "atom_indices": [1, 2],
                    "trivial_name": {
                        "name": "hydroxyl",
                        "smarts": "[OX2H]",
                        "group": "functional_group",
                        ...
                    }
                }
            ]
        }
    """
    print(f"{_COLOR_GREEN}[MCP TOOL] annotate_molecule called with smiles={smiles[:50]}{_COLOR_RESET}", flush=True)
    try:
        from rdkit import Chem
        
        # Parse the SMILES string
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"{_COLOR_YELLOW}[MCP TOOL] annotate_molecule: Invalid SMILES{_COLOR_RESET}", flush=True)
            return {
                "error": f"Invalid SMILES string: {smiles}",
                "smiles": smiles
            }
        
        # Annotate the molecule
        result = smart_chemist.mol_to_annotation_json(mol)
        print(f"{_COLOR_GREEN}[MCP TOOL] annotate_molecule: Found {len(result.get('matches', []))} functional groups{_COLOR_RESET}", flush=True)
        return result
    except Exception as e:
        print(f"{_COLOR_YELLOW}[MCP TOOL] annotate_molecule error: {e}{_COLOR_RESET}", flush=True)
        return {
            "error": str(e),
            "smiles": smiles
        }

@mcp.tool()
def get_smiles_from_name(chemical_name: str) -> Dict[str, Any]:
    """Get SMILES string from name"""
    print(f"{_COLOR_GREEN}[MCP TOOL] get_smiles_from_name called with chemical_name={chemical_name}{_COLOR_RESET}", flush=True)
    
    get_smile_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{chemical_name}/property/CanonicalSMILES/JSON"
    response = requests.get(get_smile_url)
    response_json = response.json()
    canonical_smiles = response_json["PropertyTable"]["Properties"][0].get("CanonicalSMILES", None)
    
    if canonical_smiles is None:
        canonical_smiles = response_json["PropertyTable"]["Properties"][0].get("ConnectivitySMILES", None)
    
    if canonical_smiles is None:
        print(f"{_COLOR_YELLOW}[MCP TOOL] get_smiles_from_name: No SMILES found{_COLOR_RESET}", flush=True)
        return {
            "error": f"No SMILES found for {chemical_name}",
        }
    
    print(f"{_COLOR_GREEN}[MCP TOOL] get_smiles_from_name: Found SMILES={canonical_smiles[:50]}{_COLOR_RESET}", flush=True)
    return {
        "smiles": canonical_smiles
    }


@mcp.tool()
def convert_identifier_to_smiles(identifier: str) -> Dict[str, Any]:
    """
    Convert various chemical identifiers to SMILES strings.
    
    This tool accepts multiple types of chemical identifiers and converts them
    to SMILES notation. It supports:
    - SMILES strings (returns as-is)
    - ChEMBL IDs (e.g., "chembl:CHEMBL50894")
    - ChEBI IDs (e.g., "chebi:138488")
    - PubChem CIDs (e.g., "pubchem:5005498")
    
    Args:
        identifier: Chemical identifier in one of the supported formats
    
    Returns:
        Dictionary containing:
        - identifier: The original input identifier
        - smiles_list: List of SMILES strings (may contain multiple for patterns)
        - count: Number of SMILES strings returned
    
    Example:
        >>> convert_identifier_to_smiles("chembl:CHEMBL25")
        {
            "identifier": "chembl:CHEMBL25",
            "smiles_list": ["CC(=O)Oc1ccccc1C(=O)O"],
            "count": 1
        }
    """
    print(f"{_COLOR_GREEN}[MCP TOOL] convert_identifier_to_smiles called with identifier={identifier}{_COLOR_RESET}", flush=True)
    try:
        smiles_list = convert_string_input_to_smiles(identifier)
        print(f"{_COLOR_GREEN}[MCP TOOL] convert_identifier_to_smiles: Converted to {len(smiles_list)} SMILES{_COLOR_RESET}", flush=True)
        return {
            "identifier": identifier,
            "smiles_list": smiles_list,
            "count": len(smiles_list)
        }
    except Exception as e:
        print(f"{_COLOR_YELLOW}[MCP TOOL] convert_identifier_to_smiles error: {e}{_COLOR_RESET}", flush=True)
        return {
            "error": str(e),
            "identifier": identifier,
            "smiles_list": [],
            "count": 0
        }

@mcp.tool()
def check_molecule_patent(smiles: str) -> Dict[str, Any]:
    """
    Check if molecule appears in patent literature (SureChEMBL).
    Quick FTO check - NOT legal advice, consult patent attorney.
    
    Args:
        smiles: SMILES string of molecule
        
    Returns:
        Patent status with is_patented, fto_status, confidence, patents list
    """
    print(f"{_COLOR_GREEN}[MCP TOOL] check_molecule_patent called with smiles={smiles[:50]}{_COLOR_RESET}", flush=True)
    try:
        result = patent_tool.check_patent(smiles)
        print(f"{_COLOR_GREEN}[MCP TOOL] check_molecule_patent: is_patented={result.get('is_patented')}{_COLOR_RESET}", flush=True)
        return result
    except ValueError as e:
        print(f"{_COLOR_YELLOW}[MCP TOOL] check_molecule_patent error: {e}{_COLOR_RESET}", flush=True)
        return {"error": str(e), "smiles": smiles, "is_patented": None}
    except Exception as e:
        print(f"{_COLOR_YELLOW}[MCP TOOL] check_molecule_patent error: {e}{_COLOR_RESET}", flush=True)
        return {"error": f"Patent check failed: {str(e)}", "smiles": smiles}


@mcp.tool()
def scan_for_bioisosteres(
    query_smiles: str,
    top_k: int = 20,
    min_similarity: float = 0.3,
    use_bio_filters: bool = True,
    logp_tolerance: float = 1.5,
    tpsa_tolerance: float = 25.0,
    h_donor_tolerance: int = 1,
    h_acceptor_tolerance: int = 2
) -> Dict[str, Any]:
    """
    Scan for bio-isosteric ring replacements for a given ring structure.
    
    This tool searches through a comprehensive database of ~34,000 ring systems
    to find bio-isosteric replacements - rings that maintain similar biological
    activity while potentially offering improved properties. The scan uses:
    - Structural similarity (Morgan fingerprints)
    - 2D pharmacophore similarity (functional group patterns)
    - Ring topology matching (ring sizes and fusion patterns)
    - Physicochemical property filters (logP, TPSA, H-bonding, aromaticity)
    
    Results are searched hierarchically from three sources (highest priority first):
    1. Ertl dataset - well-characterized medicinal chemistry rings
    2. ChemSpace collections - focused bioisostere libraries
    3. Novel ring clusters - diverse ring systems from literature
    
    Args:
        query_smiles: SMILES string of the query ring structure
                     Examples: "c1ccccc1" (benzene), "C1CCNCC1" (piperidine)
        top_k: Number of top matches to return (default: 20, range: 1-100)
        min_similarity: Minimum Tanimoto similarity threshold (default: 0.3, range: 0.0-1.0)
        use_bio_filters: Apply physicochemical property filters (default: True)
        logp_tolerance: Maximum logP difference for bio-isosteres (default: 1.5)
        tpsa_tolerance: Maximum TPSA difference in Ų (default: 25.0)
        h_donor_tolerance: Maximum H-donor count difference (default: 1)
        h_acceptor_tolerance: Maximum H-acceptor count difference (default: 2)
    
    Returns:
        Dictionary containing:
        - query_smiles: The input query SMILES
        - num_results: Number of bio-isosteres found
        - results: List of bio-isosteric replacements, each containing:
            - source: Data source ("ertl", "chemspace", or "clusters")
            - centroid_smiles: SMILES of the bio-isostere
            - similarity: Structural similarity score (0-1)
            - bio_isostere_score: Combined bio-isostere score (0-1)
            - pharmacophore_similarity: 2D pharmacophore similarity (0-1)
            - topology_similarity: Ring topology similarity (0-1)
            - descriptors: Physicochemical properties (logP, TPSA, etc.)
            - delta_properties: Property differences from query
            - cluster_id: (for cluster results) ID of the cluster
            - num_members: (for cluster results) Number of members in cluster
            - example_smiles: (for cluster results) Example SMILES from cluster
    
    Raises:
        ValueError: If the SMILES string is invalid or scan fails
    
    Example:
        >>> scan_for_bioisosteres("c1ccccc1", top_k=10, min_similarity=0.2, use_bio_filters=True)
        {
            "query_smiles": "c1ccccc1",
            "num_results": 5,
            "results": [
                {
                    "source": "ertl",
                    "centroid_smiles": "c1ccncc1",
                    "similarity": 0.85,
                    "bio_isostere_score": 0.82,
                    "pharmacophore_similarity": 0.78,
                    "topology_similarity": 0.95,
                    "descriptors": {"logp": 0.65, "tpsa": 12.9, ...},
                    "delta_properties": {"delta_logp": 1.24, ...}
                },
                ...
            ]
        }
    """
    print(f"{_COLOR_GREEN}[MCP TOOL] scan_for_bioisosteres called with query_smiles={query_smiles[:50]}, top_k={top_k}, min_similarity={min_similarity}{_COLOR_RESET}", flush=True)
    try:
        # Call the scanner
        results = ring_scanner.scan(
            query_smiles=query_smiles,
            top_k=top_k,
            min_similarity=min_similarity,
            use_bio_filters=use_bio_filters,
            logp_tolerance=logp_tolerance,
            tpsa_tolerance=tpsa_tolerance,
            h_donor_tolerance=h_donor_tolerance,
            h_acceptor_tolerance=h_acceptor_tolerance
        )
        
        print(f"{_COLOR_GREEN}[MCP TOOL] scan_for_bioisosteres: Found {len(results)} bioisosteres{_COLOR_RESET}", flush=True)
        return {
            "query_smiles": query_smiles,
            "num_results": len(results),
            "results": results
        }
    except Exception as e:
        print(f"{_COLOR_YELLOW}[MCP TOOL] scan_for_bioisosteres error: {e}{_COLOR_RESET}", flush=True)
        return {
            "error": str(e),
            "query_smiles": query_smiles,
            "num_results": 0,
            "results": []
        }


def main():
    mcp.run()

if __name__ == "__main__":
    ivermectin_smiles = "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"
    output = annotate_molecule(ivermectin_smiles)
    print(output)
