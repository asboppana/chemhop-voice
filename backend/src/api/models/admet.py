"""
Request and response models for ADMET prediction endpoints.
"""
from typing import Dict, Any
from pydantic import BaseModel, Field


class ADMETRequest(BaseModel):
    """Request model for ADMET predictions.
    
    - smiles: SMILES string representation of the molecule
    """
    smiles: str = Field(
        ...,
        description="SMILES string representation of the molecule",
        min_length=1,
        examples=["O(c1ccc(cc1)CCOC)CC(O)CNC(C)C", "CCO", "c1ccccc1"]
    )


class ADMETResponse(BaseModel):
    """Response model for ADMET predictions.
    
    Contains the predicted ADMET properties for the molecule.
    The predictions dictionary contains various pharmacokinetic and 
    toxicity properties predicted by the ADMET AI model.
    """
    smiles: str = Field(..., description="Input SMILES string")
    predictions: Dict[str, Any] = Field(
        ..., 
        description="Dictionary of predicted ADMET properties and their values"
    )

