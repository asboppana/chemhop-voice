"""
ADMET prediction endpoints.

Provides functionality for predicting ADMET (Absorption, Distribution, 
Metabolism, Excretion, and Toxicity) properties of molecules using SMILES strings.
"""
from fastapi import APIRouter, Depends, status, HTTPException

from src.api.models import ErrorResponse
from src.api.models.admet import ADMETRequest, ADMETResponse
from src.controllers.admet_controller import ADMETController


# ============================================================================
# Dependency Injection
# ============================================================================


def get_admet_controller() -> ADMETController:
    """Dependency injection for ADMETController."""
    return ADMETController()


# ============================================================================
# Router
# ============================================================================

router = APIRouter()


# ============================================================================
# Endpoints
# ============================================================================


@router.post(
    "/admet/predict",
    status_code=status.HTTP_200_OK,
    response_model=ADMETResponse,
    responses={
        400: {"model": ErrorResponse, "description": "Invalid request"},
        422: {"model": ErrorResponse, "description": "Invalid SMILES string or prediction failed"},
        500: {"model": ErrorResponse, "description": "Internal server error"},
    },
)
async def predict_admet(
    request: ADMETRequest,
    controller: ADMETController = Depends(get_admet_controller),
) -> ADMETResponse:
    """
    Predict ADMET properties for a molecule from its SMILES string.
    
    Takes a SMILES string representation of a molecule and returns:
    - Original SMILES string
    - Dictionary of predicted ADMET properties including:
      - Absorption properties
      - Distribution properties
      - Metabolism properties
      - Excretion properties
      - Toxicity properties
    
    The endpoint uses the ADMET AI model to predict various pharmacokinetic
    and toxicity properties that are important for drug development.
    """
    try:
        result = controller.predict_admet(request.smiles)
        return result
    except ValueError as e:
        raise HTTPException(
            status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
            detail=str(e)
        )
    except Exception as e:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Error predicting ADMET properties: {str(e)}"
        )

