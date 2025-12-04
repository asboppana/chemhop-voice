"""
Molecule analysis endpoints.

Provides functionality for analyzing molecules using SMILES strings
and returning annotated chemical structure information.
"""
from fastapi import APIRouter, Depends, status, HTTPException
from pydantic import BaseModel, Field

from src.api.models import ErrorResponse
from src.api.models.molecule import (
    MoleculeRequest,
    MoleculeResponse,
    HighlightRequest,
    HighlightResponse,
    BioisostereScanRequest,
    BioisostereScanResponse,
    ClusterInfoResponse,
    ExtractSubstructureRequest,
    ExtractSubstructureResponse,
    GenerateSvgRequest,
    GenerateSvgResponse,
    ReplaceSubstructureRequest,
    ReplaceSubstructureResponse,
    LLMReplaceSubstructureRequest,
    LLMReplaceSubstructureResponse
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


@router.post(
    "/molecule/extract-substructure",
    status_code=status.HTTP_200_OK,
    response_model=ExtractSubstructureResponse,
    responses={
        400: {"model": ErrorResponse, "description": "Invalid request"},
        422: {"model": ErrorResponse, "description": "Invalid SMILES or atom indices"},
        500: {"model": ErrorResponse, "description": "Internal server error"},
    },
)
async def extract_substructure(
    request: ExtractSubstructureRequest,
    controller: MoleculeController = Depends(get_molecule_controller),
) -> ExtractSubstructureResponse:
    """
    Extract a substructure SMILES from a parent molecule.
    
    Takes a parent molecule SMILES and a list of atom indices, then extracts
    the substructure defined by those atoms as a separate SMILES string.
    
    Parameters:
    - **smiles**: SMILES string of the parent molecule
    - **atom_indices**: List of atom indices that define the substructure
    
    Returns:
    - Substructure SMILES
    - Parent SMILES
    """
    try:
        substructure_smiles = controller.extract_substructure_smiles(
            request.smiles,
            request.atom_indices
        )
        return ExtractSubstructureResponse(
            substructure_smiles=substructure_smiles,
            parent_smiles=request.smiles
        )
    except ValueError as e:
        raise HTTPException(
            status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
            detail=str(e)
        )
    except Exception as e:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Error extracting substructure: {str(e)}"
        )


@router.post(
    "/molecule/generate-svg",
    status_code=status.HTTP_200_OK,
    response_model=GenerateSvgResponse,
    responses={
        400: {"model": ErrorResponse, "description": "Invalid request"},
        422: {"model": ErrorResponse, "description": "Invalid SMILES string"},
        500: {"model": ErrorResponse, "description": "Internal server error"},
    },
)
async def generate_svg(
    request: GenerateSvgRequest,
    controller: MoleculeController = Depends(get_molecule_controller),
) -> GenerateSvgResponse:
    """
    Generate an SVG visualization from a SMILES string.
    
    Takes a SMILES string and generates an SVG representation of the molecule.
    
    Parameters:
    - **smiles**: SMILES string of the molecule
    - **width**: Width of the SVG in pixels (default: 200)
    - **height**: Height of the SVG in pixels (default: 200)
    
    Returns:
    - SVG string
    - Original SMILES
    """
    try:
        svg = controller.generate_molecule_svg(
            request.smiles,
            request.width,
            request.height
        )
        return GenerateSvgResponse(
            svg=svg,
            smiles=request.smiles
        )
    except ValueError as e:
        raise HTTPException(
            status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
            detail=str(e)
        )
    except Exception as e:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Error generating SVG: {str(e)}"
        )


@router.post(
    "/molecule/replace-substructure",
    status_code=status.HTTP_200_OK,
    response_model=ReplaceSubstructureResponse,
    responses={
        400: {"model": ErrorResponse, "description": "Invalid request"},
        422: {"model": ErrorResponse, "description": "Invalid SMILES or replacement failed"},
        500: {"model": ErrorResponse, "description": "Internal server error"},
    },
)
async def replace_substructure(
    request: ReplaceSubstructureRequest,
    controller: MoleculeController = Depends(get_molecule_controller),
) -> ReplaceSubstructureResponse:
    """
    Replace a substructure in a molecule with a bio-isostere replacement.
    
    Attempts to optimally align the replacement pattern with the source pattern,
    trying 3 different rotations (0°, 120°, 240°) to minimize distance between
    connection points on the perimeter.
    
    Parameters:
    - **parent_smiles**: SMILES of the parent molecule
    - **source_pattern_smiles**: SMILES of the pattern to replace
    - **replacement_smiles**: SMILES of the replacement pattern
    - **atom_indices**: Atom indices of the source pattern in the parent molecule
    
    Returns:
    - Original parent SMILES
    - List of generated molecule SMILES (up to 3)
    - Number of successfully generated molecules
    """
    try:
        replacement_smiles_list = controller.replace_substructure_with_alignment(
            request.parent_smiles,
            request.source_pattern_smiles,
            request.replacement_smiles,
            request.atom_indices
        )
        
        return ReplaceSubstructureResponse(
            original_smiles=request.parent_smiles,
            replacement_smiles_list=replacement_smiles_list,
            num_generated=len(replacement_smiles_list)
        )
    except ValueError as e:
        raise HTTPException(
            status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
            detail=str(e)
        )
    except Exception as e:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Error replacing substructure: {str(e)}"
        )


@router.post(
    "/molecule/llm-replace-substructure",
    response_model=LLMReplaceSubstructureResponse,
    responses={
        status.HTTP_200_OK: {
            "model": LLMReplaceSubstructureResponse,
            "description": "Successful LLM-based substructure replacement"
        },
        status.HTTP_422_UNPROCESSABLE_ENTITY: {
            "model": ErrorResponse,
            "description": "Invalid SMILES strings"
        },
        status.HTTP_500_INTERNAL_SERVER_ERROR: {
            "model": ErrorResponse,
            "description": "Internal server error during replacement"
        }
    },
    summary="Replace substructure using LLM",
    description="Use GPT to intelligently replace a fragment in a molecule with another fragment"
)
async def llm_replace_substructure(
    request: LLMReplaceSubstructureRequest,
    controller: MoleculeController = Depends(get_molecule_controller)
) -> LLMReplaceSubstructureResponse:
    """
    Replace a substructure in a molecule using GPT.
    
    This endpoint takes:
    - **original_smiles**: SMILES of the original molecule
    - **source_fragment_smiles**: SMILES of the fragment to be replaced
    - **replacement_fragment_smiles**: SMILES of the replacement fragment
    
    Returns:
    - The resulting molecule SMILES after replacement
    - Explanation from the LLM
    - Success status
    """
    try:
        result = controller.replace_substructure_with_llm(
            request.original_smiles,
            request.source_fragment_smiles,
            request.replacement_fragment_smiles
        )
        
        return LLMReplaceSubstructureResponse(
            original_smiles=request.original_smiles,
            source_fragment_smiles=request.source_fragment_smiles,
            replacement_fragment_smiles=request.replacement_fragment_smiles,
            result_smiles=result["result_smiles"],
            explanation=result.get("explanation"),
            success=result["success"]
        )
    except ValueError as e:
        raise HTTPException(
            status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
            detail=str(e)
        )
    except Exception as e:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Error in LLM replacement: {str(e)}"
        )


# ============================================================================
# ADMET Prediction Endpoint
# ============================================================================

class ADMETPredictionRequest(BaseModel):
    """Request model for ADMET prediction."""
    smiles: str = Field(
        ...,
        description="SMILES string representation of the molecule",
        min_length=1
    )


class ADMETPredictionResponse(BaseModel):
    """Response model for ADMET prediction."""
    smiles: str = Field(..., description="Input SMILES string")
    predictions: dict = Field(default_factory=dict, description="ADMET property predictions")
    error: str | None = Field(None, description="Error message if prediction failed")


@router.post(
    "/molecule/predict-admet",
    response_model=ADMETPredictionResponse,
    responses={
        status.HTTP_200_OK: {
            "model": ADMETPredictionResponse,
            "description": "Successful ADMET prediction"
        },
        status.HTTP_422_UNPROCESSABLE_ENTITY: {
            "model": ErrorResponse,
            "description": "Invalid SMILES string"
        },
        status.HTTP_500_INTERNAL_SERVER_ERROR: {
            "model": ErrorResponse,
            "description": "Internal server error during prediction"
        }
    },
    summary="Predict ADMET properties",
    description="Predict ADMET (Absorption, Distribution, Metabolism, Excretion, Toxicity) properties for a molecule"
)
async def predict_admet(
    request: ADMETPredictionRequest
) -> ADMETPredictionResponse:
    """
    Predict ADMET properties for a molecule.
    
    This endpoint calls the ADMET AI service to predict pharmacokinetic
    and toxicity properties including:
    - Absorption: Caco-2, Solubility, HIA, Pgp
    - Distribution: BBB, PPB, VDss
    - Metabolism: CYP inhibition
    - Excretion: Half-life, Clearance
    - Toxicity: hERG, AMES, DILI
    """
    import requests
    import os
    
    try:
        # Call the ADMET service
        admet_service_url = os.getenv("ADMET_SERVICE_URL", "http://localhost:8002")
        
        # Try to call the MCP tool endpoint
        response = requests.post(
            f"{admet_service_url}/call-tool",
            json={
                "name": "predict_admet_properties",
                "arguments": {"smiles": request.smiles}
            },
            timeout=30
        )
        
        if response.status_code == 200:
            result = response.json()
            # MCP response format has content array
            if "content" in result and len(result["content"]) > 0:
                content = result["content"][0]
                if "text" in content:
                    import json
                    predictions = json.loads(content["text"])
                    return ADMETPredictionResponse(
                        smiles=request.smiles,
                        predictions=predictions.get("predictions", predictions),
                        error=predictions.get("error")
                    )
            
            return ADMETPredictionResponse(
                smiles=request.smiles,
                predictions=result.get("predictions", {}),
                error=result.get("error")
            )
        else:
            return ADMETPredictionResponse(
                smiles=request.smiles,
                predictions={},
                error=f"ADMET service returned status {response.status_code}"
            )
            
    except requests.exceptions.ConnectionError:
        return ADMETPredictionResponse(
            smiles=request.smiles,
            predictions={},
            error="ADMET service not available"
        )
    except Exception as e:
        return ADMETPredictionResponse(
            smiles=request.smiles,
            predictions={},
            error=str(e)
        )
