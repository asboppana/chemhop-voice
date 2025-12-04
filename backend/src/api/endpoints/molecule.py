"""
Molecule analysis endpoints.

Provides functionality for analyzing molecules using SMILES strings
and returning annotated chemical structure information.
"""
from fastapi import APIRouter, Depends, status, HTTPException

from src.api.models import ErrorResponse
from src.api.models.molecule import (
    MoleculeRequest,
    MoleculeResponse,
    HighlightRequest,
    HighlightResponse,
    BioisostereScanRequest,
    BioisostereScanResponse,
    ClusterInfoResponse
)
from src.controllers.molecule_controller import MoleculeController
from src.controllers.bioisostere_controller import BioisostereController

# ============================================================================
# Dependency Injection
# ============================================================================


def get_molecule_controller() -> MoleculeController:
    """Dependency injection for MoleculeController."""
    return MoleculeController()


def get_bioisostere_controller() -> BioisostereController:
    """Dependency injection for BioisostereController."""
    return BioisostereController()


# ============================================================================
# Router
# ============================================================================

router = APIRouter()


# ============================================================================
# Endpoints
# ============================================================================


@router.post(
    "/molecule/analyze",
    status_code=status.HTTP_200_OK,
    response_model=MoleculeResponse,
    responses={
        400: {"model": ErrorResponse, "description": "Invalid request"},
        422: {"model": ErrorResponse, "description": "Invalid SMILES string"},
        500: {"model": ErrorResponse, "description": "Internal server error"},
    },
)
async def analyze_molecule(
    request: MoleculeRequest,
    controller: MoleculeController = Depends(get_molecule_controller),
) -> MoleculeResponse:
    """
    Analyze a molecule from its SMILES string.
    
    Takes a SMILES string representation of a molecule and returns:
    - Molecule name (if available)
    - SVG visualization of the molecule structure
    - List of matched chemical patterns/functional groups
    - Original SMILES string
    
    The endpoint uses the SmartChemist tool to identify and annotate
    various chemical substructures and functional groups present in the molecule.
    """
    try:
        result = controller.analyze_molecule(request.smiles)
        return result
    except ValueError as e:
        raise HTTPException(
            status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
            detail=f"Invalid SMILES string: {str(e)}"
        )
    except Exception as e:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Error analyzing molecule: {str(e)}"
        )


@router.post(
    "/molecule/highlight",
    status_code=status.HTTP_200_OK,
    response_model=HighlightResponse,
    responses={
        400: {"model": ErrorResponse, "description": "Invalid request"},
        422: {"model": ErrorResponse, "description": "Invalid SMILES string"},
        500: {"model": ErrorResponse, "description": "Internal server error"},
    },
)
async def highlight_molecule(
    request: HighlightRequest,
    controller: MoleculeController = Depends(get_molecule_controller),
) -> HighlightResponse:
    """
    Generate an SVG of a molecule with specific atoms highlighted.
    
    Takes a SMILES string and a list of atom indices to highlight,
    then returns an SVG visualization with those atoms highlighted.
    """
    try:
        result = controller.highlight_molecule(request.smiles, request.atom_indices)
        return result
    except ValueError as e:
        raise HTTPException(
            status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
            detail=f"Invalid SMILES string: {str(e)}"
        )
    except Exception as e:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Error highlighting molecule: {str(e)}"
        )


@router.post(
    "/molecule/bioisostere/scan",
    status_code=status.HTTP_200_OK,
    response_model=BioisostereScanResponse,
    responses={
        400: {"model": ErrorResponse, "description": "Invalid request"},
        422: {"model": ErrorResponse, "description": "Invalid SMILES string"},
        500: {"model": ErrorResponse, "description": "Internal server error"},
    },
)
async def scan_bioisosteres(
    request: BioisostereScanRequest,
    controller: BioisostereController = Depends(get_bioisostere_controller),
) -> BioisostereScanResponse:
    """
    Scan for bio-isosteric ring replacements with adaptive scoring.
    
    **Hierarchical Data Sources:**
    - ertl_data: 189 curated bio-isosteres (highest priority)
    - Chemspace: 10,474 medicinal chemistry scaffolds
    - Clustered DB: 249,281 centroids representing ~4M rings
    
    **Adaptive Scoring (Pure Scoring - No Filters):**
    - Small molecules (<10 atoms): Pharmacophore-driven (20% Morgan, 50% Pharm, 30% Topo)
    - Large molecules (≥10 atoms): Balanced approach (40% Morgan, 40% Pharm, 20% Topo)
    - All candidates ranked by composite score
    - Property differences provided for reference
    
    **Parameters:**
    - **ring_smiles**: SMILES string of the query ring structure
    - **top_k**: Number of top matches to return (1-100, default: 20)
    - **min_similarity**: Minimum Morgan FP similarity threshold (0.0-1.0, default: 0.3)
    
    **Returns:**
    - Query SMILES
    - Number of matches found
    - List of bio-isosteres with:
      - **bio_isostere_score**: Composite adaptive score (0-1)
      - **similarity**: Morgan fingerprint Tanimoto similarity
      - **pharmacophore_similarity**: Gobbi 2D pharmacophore score
      - **topology_similarity**: Ring architecture score
      - **source**: Data source (ertl, chemspace, or clusters)
      - **descriptors**: Full physicochemical profile
      - **delta_properties**: Property differences from query (logP, TPSA, H-bonding)
    """
    try:
        result = controller.scan_ring(
            ring_smiles=request.ring_smiles,
            top_k=request.top_k,
            min_similarity=request.min_similarity
        )
        return result
    except ValueError as e:
        raise HTTPException(
            status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
            detail=str(e)
        )
    except Exception as e:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Error scanning for bio-isosteres: {str(e)}"
        )


@router.get(
    "/molecule/bioisostere/cluster/{cluster_id}",
    status_code=status.HTTP_200_OK,
    response_model=ClusterInfoResponse,
    responses={
        404: {"model": ErrorResponse, "description": "Cluster not found"},
        500: {"model": ErrorResponse, "description": "Internal server error"},
    },
)
async def get_cluster_info(
    cluster_id: int,
    controller: BioisostereController = Depends(get_bioisostere_controller),
) -> ClusterInfoResponse:
    """
    Get detailed information about a specific cluster.
    
    Retrieves all member SMILES and statistics for a cluster identified
    during bio-isostere scanning.
    
    Parameters:
    - **cluster_id**: ID of the cluster to retrieve information about
    
    Returns:
    - Cluster ID
    - Number of members in the cluster
    - List of all SMILES strings in the cluster
    """
    try:
        result = controller.get_cluster_details(cluster_id)
        return result
    except ValueError as e:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=str(e)
        )
    except Exception as e:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Error retrieving cluster info: {str(e)}"
        )
