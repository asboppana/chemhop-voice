"""
ADMET AI MCP Server.

This server exposes ADMET property prediction via the Model Context Protocol (MCP).
It provides access to ADMET AI model for predicting molecular properties.

Run with: python main.py
"""
from typing import Any, Dict

from fastmcp import FastMCP

from src.admet_tool import ADMETAITool

# Initialize FastMCP server
mcp = FastMCP("admet_mcp")

# Initialize ADMET tool instance
admet_tool = ADMETAITool()


@mcp.tool()
def predict_admet_properties(smiles: str) -> Dict[str, Any]:
    """
    Predict ADMET (Absorption, Distribution, Metabolism, Excretion, and Toxicity) 
    properties for a molecule using its SMILES string representation.
    
    This tool uses the ADMET AI model to predict various pharmacokinetic and 
    toxicity properties that are critical for drug development, including:
    - Absorption properties (e.g., Caco-2 permeability, solubility)
    - Distribution properties (e.g., plasma protein binding, volume of distribution)
    - Metabolism properties (e.g., CYP enzyme inhibition)
    - Excretion properties (e.g., clearance, half-life)
    - Toxicity properties (e.g., hERG inhibition, hepatotoxicity)
    
    Args:
        smiles: SMILES string representation of the molecule
                Examples: "CCO" (ethanol), "c1ccccc1" (benzene), 
                         "CC(=O)Oc1ccccc1C(=O)O" (aspirin)
    
    Returns:
        Dictionary containing:
        - smiles: The input SMILES string
        - predictions: Dictionary of predicted ADMET properties with property 
                      names as keys and predicted values/scores as values
    
    Raises:
        ValueError: If the SMILES string is invalid or prediction fails
    
    Example:
        >>> predict_admet_properties("CCO")
        {
            "smiles": "CCO",
            "predictions": {
                "Caco2_Wang": 0.234,
                "Solubility_AqSolDB": -0.234,
                "HIA_Hou": 0.95,
                ...
            }
        }
    """
    try:
        predictions = admet_tool.predict(smiles)
        return {
            "smiles": smiles,
            "predictions": predictions
        }
    except Exception as e:
        return {
            "error": str(e),
            "smiles": smiles
        }


def run():
    """Run the MCP server."""
    mcp.run()

