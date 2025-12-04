"""
Controller for ADMET prediction operations.

Handles the business logic for predicting ADMET properties of molecules
using SMILES strings and interfacing with the ADMET AI tool.
"""
from src.services.mcp.admet_tools.admet_ai_tool import ADMETAITool
from src.api.models.admet import ADMETResponse


class ADMETController:
    """Controller for ADMET prediction operations."""
    
    def __init__(self):
        """Initialize the ADMETController with an ADMET AI tool instance."""
        self.admet_tool = ADMETAITool()
    
    def predict_admet(self, smiles: str) -> ADMETResponse:
        """
        Predict ADMET properties for a molecule from its SMILES string.
        
        Args:
            smiles: SMILES string representation of the molecule
            
        Returns:
            ADMETResponse containing the predicted ADMET properties
            
        Raises:
            ValueError: If the SMILES string is invalid or prediction fails
        """
        # Get predictions from ADMET AI tool
        predictions = self.admet_tool.predict(smiles)
        
        # Create and return the response
        return ADMETResponse(
            smiles=smiles,
            predictions=predictions
        )

