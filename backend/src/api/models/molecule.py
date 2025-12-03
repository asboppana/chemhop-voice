"""
Request and response models for molecule analysis endpoints.
"""
from typing import List, Dict, Any, Optional
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


class BioisostereScanRequest(BaseModel):
    """Request model for bio-isostere scanning.
    
    Uses pure adaptive scoring to rank bio-isosteres:
    - Small molecules (<10 atoms): Pharmacophore-driven (20% Morgan, 50% Pharm, 30% Topo)
    - Large molecules (≥10 atoms): Balanced (40% Morgan, 40% Pharm, 20% Topo)
    
    - ring_smiles: SMILES string of the query ring to find replacements for
    - top_k: Number of top matches to return
    - min_similarity: Minimum Morgan fingerprint similarity threshold
    """
    ring_smiles: str = Field(
        ...,
        description="SMILES string of the query ring",
        min_length=1,
        examples=["c1ccccc1", "c1ccc2ccccc2c1", "C1CCCCC1"]
    )
    top_k: int = Field(
        20,
        description="Number of top matches to return",
        ge=1,
        le=100
    )
    min_similarity: float = Field(
        0.3,
        description="Minimum Morgan fingerprint similarity threshold (0.0-1.0)",
        ge=0.0,
        le=1.0
    )


class BioisostereMatch(BaseModel):
    """A matched bio-isosteric ring with adaptive scoring."""
    source: str = Field(..., description="Data source (ertl, chemspace, or clusters)")
    similarity: float = Field(..., description="Morgan fingerprint Tanimoto similarity")
    centroid_smiles: str = Field(..., description="Representative SMILES")
    bio_isostere_score: float = Field(..., description="Composite bio-isostere score (0-1)")
    pharmacophore_similarity: float = Field(..., description="Gobbi 2D pharmacophore similarity")
    topology_similarity: float = Field(..., description="Ring topology similarity")
    cluster_id: Optional[int] = Field(None, description="Cluster ID (only for cluster source)")
    num_members: int = Field(0, description="Number of rings in cluster (0 for non-cluster)")
    example_smiles: List[str] = Field(default_factory=list, description="Example SMILES (for clusters)")
    descriptors: Dict[str, Any] = Field(default_factory=dict, description="Physicochemical descriptors")
    delta_properties: Dict[str, float] = Field(default_factory=dict, description="Property differences from query")


class BioisostereScanResponse(BaseModel):
    """Response model for bio-isostere scanning.
    
    Contains matches of similar ring systems that could serve as bio-isosteric replacements.
    """
    query_smiles: str = Field(..., description="Original query ring SMILES")
    num_matches: int = Field(..., description="Number of matches found")
    matches: List[BioisostereMatch] = Field(..., description="List of matched bio-isosteric clusters")


class ClusterInfoRequest(BaseModel):
    """Request model for getting detailed cluster information."""
    cluster_id: int = Field(
        ...,
        description="ID of the cluster to get information about",
        ge=0
    )


class ClusterInfoResponse(BaseModel):
    """Response model for cluster information."""
    cluster_id: int = Field(..., description="ID of the cluster")
    num_members: int = Field(..., description="Number of rings in this cluster")
    member_smiles: List[str] = Field(..., description="All SMILES in the cluster")
