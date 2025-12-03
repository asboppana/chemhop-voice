"""
Molecule analysis endpoints.

Provides functionality for analyzing molecules using SMILES strings
and returning annotated chemical structure information.
"""
from fastapi import APIRouter, Depends, status, HTTPException

from src.api.models import ErrorResponse
from src.api.models.molecule import MoleculeRequest, MoleculeResponse
from src.controllers.molecule_controller import MoleculeController

# ============================================================================
# Dependency Injection
# ============================================================================


def get_molecule_controller() -> MoleculeController:
    """Dependency injection for MoleculeController."""
    return MoleculeController()


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

