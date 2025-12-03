"""
Controller for bio-isostere scanning operations.

Handles the business logic for scanning ring structures with adaptive scoring
to find bio-isosteric replacements.
"""
import logging
from typing import Optional

from src.services.mcp.ring_scan import get_scanner
from src.api.models.molecule import (
    BioisostereScanResponse,
    BioisostereMatch,
    ClusterInfoResponse
)

logger = logging.getLogger(__name__)


class BioisostereController:
    """Controller for bio-isostere scanning operations."""
    
    def __init__(self):
        """Initialize the BioisostereController with a scanner instance."""
        self.scanner = get_scanner()
    
    def scan_ring(
        self,
        ring_smiles: str,
        top_k: int = 20,
        min_similarity: float = 0.3
    ) -> BioisostereScanResponse:
        """
        Scan for bio-isosteric replacements of a ring structure.
        
        Uses hierarchical scanning with adaptive weighting:
        - Small molecules (<10 atoms): Pharmacophore-driven (20% Morgan, 50% Pharm, 30% Topo)
        - Large molecules (≥10 atoms): Balanced approach (40% Morgan, 40% Pharm, 20% Topo)
        
        Pure scoring approach - no hard filters, all candidates ranked by composite score.
        
        Args:
            ring_smiles: SMILES string of the query ring
            top_k: Number of top matches to return
            min_similarity: Minimum Morgan FP similarity threshold
            
        Returns:
            BioisostereScanResponse containing the matches
            
        Raises:
            ValueError: If the SMILES string is invalid
        """
        logger.info(
            f"Scanning ring: {ring_smiles} "
            f"(top_k={top_k}, min_similarity={min_similarity})"
        )
        
        try:
            # Perform the scan with adaptive scoring
            results = self.scanner.scan(
                query_smiles=ring_smiles,
                top_k=top_k,
                min_similarity=min_similarity,
                use_bio_filters=True  # Always calculate descriptors for scoring
            )
            
            # Convert to response models
            matches = [
                BioisostereMatch(
                    source=r.get("source", "unknown"),
                    similarity=r["similarity"],
                    centroid_smiles=r["centroid_smiles"],
                    bio_isostere_score=r["bio_isostere_score"],
                    pharmacophore_similarity=r["pharmacophore_similarity"],
                    topology_similarity=r["topology_similarity"],
                    cluster_id=r.get("cluster_id"),
                    num_members=r.get("num_members", 0),
                    example_smiles=r.get("example_smiles", []),
                    descriptors=r.get("descriptors", {}),
                    delta_properties=r.get("delta_properties", {})
                )
                for r in results
            ]
            
            logger.info(
                f"Found {len(matches)} bio-isosteres for ring {ring_smiles} "
                f"(sources: {self._count_sources(matches)})"
            )
            
            return BioisostereScanResponse(
                query_smiles=ring_smiles,
                num_matches=len(matches),
                matches=matches
            )
            
        except ValueError as e:
            logger.error(f"Invalid SMILES string: {ring_smiles}")
            raise ValueError(f"Invalid SMILES string: {str(e)}")
        except Exception as e:
            logger.error(f"Error scanning ring {ring_smiles}: {str(e)}")
            raise
    
    def get_cluster_details(self, cluster_id: int) -> ClusterInfoResponse:
        """
        Get detailed information about a specific cluster.
        
        Args:
            cluster_id: ID of the cluster
            
        Returns:
            ClusterInfoResponse with cluster details
            
        Raises:
            ValueError: If the cluster ID is invalid
        """
        logger.info(f"Getting info for cluster {cluster_id}")
        
        cluster_info = self.scanner.get_cluster_info(cluster_id)
        
        if cluster_info is None:
            raise ValueError(f"Cluster ID {cluster_id} not found")
        
        return ClusterInfoResponse(
            cluster_id=cluster_info["cluster_id"],
            num_members=cluster_info["num_members"],
            member_smiles=cluster_info["member_smiles"]
        )
    
    def _count_sources(self, matches: list) -> str:
        """Count matches by source for logging."""
        sources = {}
        for match in matches:
            source = match.source
            sources[source] = sources.get(source, 0) + 1
        
        return ", ".join(f"{k}={v}" for k, v in sources.items())
