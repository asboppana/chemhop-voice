"""
Request and response models for molecule analysis endpoints.
"""
from typing import List, Dict, Any
from pydantic import BaseModel, Field


class MoleculeRequest(BaseModel):
    """Request model for molecule analysis.
    
    - smiles: SMILES string representation of the molecule
    """
    smiles: str = Field(
        ...,
        description="SMILES string representation of the molecule",
        min_length=1,
        examples=["CCO", "c1ccccc1", "CC(=O)O"]
    )


class TrivialName(BaseModel):
    """Information about a matched chemical pattern."""
    name: str = Field(..., description="Common name of the chemical pattern")
    smarts: str = Field(..., description="SMARTS pattern string")
    group: str = Field(..., description="Classification group of the pattern")
    bonds: int = Field(..., description="Number of bonds in the pattern")
    hierarchy: str = Field(..., description="Hierarchical relationship information")
    index: int = Field(..., description="Index identifier for the pattern")


class Match(BaseModel):
    """A matched chemical pattern in the molecule."""
    atom_indices: List[int] = Field(..., description="List of atom indices that match the pattern")
    trivial_name: TrivialName = Field(..., description="Information about the matched pattern")


class MoleculeResponse(BaseModel):
    """Response model for molecule analysis.
    
    Contains the annotated molecule information including visualization
    and identified chemical patterns.
    """
    name: str = Field(..., description="Name of the molecule")
    svg: str = Field(..., description="SVG representation of the molecule structure")
    matches: List[Match] = Field(..., description="List of matched chemical patterns")
    smiles: str = Field(..., description="Original SMILES string")

