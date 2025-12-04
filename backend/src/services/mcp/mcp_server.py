"""
MCP Server for Drug Discovery Tools.

This server exposes chemistry and drug discovery tools via the Model Context Protocol (MCP).
It provides access to:
- Smart Chemist: Annotate molecular structures with functional groups

Run with: python -m src.services.mcp.mcp_server
"""

import sys
from pathlib import Path
from typing import Any, Dict

from fastmcp import FastMCP

# Add the backend directory to sys.path for imports
backend_dir = Path(__file__).resolve().parent.parent.parent.parent
if str(backend_dir) not in sys.path:
    sys.path.insert(0, str(backend_dir))

# Import our tools
from src.services.mcp.smart_chemist.tools.smart_chemist import (
    SmartChemist,
    convert_string_input_to_smiles,
)
from src.services.mcp.patent_search.surechembl_tool import SureChEMBLPatentTool

# Initialize FastMCP server
mcp = FastMCP("drugdiscovery_mcp")

# Initialize tool instances
smart_chemist = SmartChemist()
patent_tool = SureChEMBLPatentTool()

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
    try:
        from rdkit import Chem
        
        # Parse the SMILES string
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {
                "error": f"Invalid SMILES string: {smiles}",
                "smiles": smiles
            }
        
        # Annotate the molecule
        result = smart_chemist.mol_to_annotation_json(mol)
        return result
    except Exception as e:
        return {
            "error": str(e),
            "smiles": smiles
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
    try:
        smiles_list = convert_string_input_to_smiles(identifier)
        return {
            "identifier": identifier,
            "smiles_list": smiles_list,
            "count": len(smiles_list)
        }
    except Exception as e:
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
    try:
        return patent_tool.check_patent(smiles)
    except ValueError as e:
        return {"error": str(e), "smiles": smiles, "is_patented": None}
    except Exception as e:
        # logger.error(f"Patent check error: {e}")
        return {"error": f"Patent check failed: {str(e)}", "smiles": smiles}


def main():
    mcp.run()

if __name__ == "__main__":
    main()