"""
ADMET AI Tool for predicting molecular properties.

This tool uses the ADMET AI model to predict various pharmacokinetic
and toxicity properties for molecules based on their SMILES representation.
"""
from typing import Dict, Any
from admet_ai import ADMETModel


class ADMETAITool:
    """
    Tool for computing ADMET (Absorption, Distribution, Metabolism, Excretion, and Toxicity) 
    predictions for molecules using the ADMET AI model.
    """
    
    def __init__(self):
        """Initialize the ADMET AI model."""
        self.model = ADMETModel()
    
    def predict(self, smiles: str) -> Dict[str, Any]:
        """
        Predict ADMET properties for a molecule.
        
        Args:
            smiles: SMILES string representation of the molecule
            
        Returns:
            Dictionary containing ADMET predictions with property names as keys
            and predicted values/scores
            
        Raises:
            ValueError: If the SMILES string is invalid or prediction fails
        """
        if not smiles or not isinstance(smiles, str):
            raise ValueError("SMILES string must be a non-empty string")
        
        try:
            # Get predictions from ADMET AI model
            preds = self.model.predict(smiles=smiles)
            
            # Convert predictions to a serializable format
            # The model returns a DataFrame-like object
            result = {}
            
            # Extract predictions - the exact format may vary based on model version
            # but typically returns a dictionary or DataFrame with property predictions
            if hasattr(preds, 'to_dict'):
                # If it's a DataFrame
                result = preds.to_dict('records')[0] if len(preds) > 0 else {}
            elif isinstance(preds, dict):
                result = preds
            else:
                # Try to convert to dict if possible
                try:
                    result = dict(preds)
                except (TypeError, ValueError):
                    result = {"predictions": str(preds)}
            
            return result
            
        except Exception as e:
            raise ValueError(f"Failed to predict ADMET properties: {str(e)}")

