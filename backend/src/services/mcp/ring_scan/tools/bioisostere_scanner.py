"""
Bio-isostere Scanner Service

Scans input ring SMILES against ~34,000 cluster representatives from
the clustered ring database to find similar ring systems that could serve
as bio-isosteric replacements.
"""
import pickle
import numpy as np
from pathlib import Path
from typing import List, Tuple, Optional
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Crippen, Lipinski, rdMolDescriptors, rdFingerprintGenerator
from rdkit.Chem.Pharm2D import Generate, Gobbi_Pharm2D
import logging

logger = logging.getLogger(__name__)


class BioisostereScanner:
    """
    Scanner for finding bio-isosteric ring replacements.
    
    Uses both structural similarity (Morgan fingerprints) and physicochemical
    properties (logP, TPSA, H-bonding, aromaticity) to identify true bio-isosteres.
    """
    
    def __init__(self):
        """Initialize the scanner by loading cluster centroids and ring data."""
        # Updated path: now we're in ring_scan/tools, so ring_data is ../ring_data
        self.ring_data_path = Path(__file__).parent.parent / "ring_data"
        self.cluster_path = self.ring_data_path / "novel_rings" / "clusters"
        self.chemspace_path = self.ring_data_path / "simple_rings"
        
        # Legacy db_path for rings.smi
        self.db_path = self.ring_data_path / "novel_rings"
        
        # Load ertl data (highest priority)
        self.ertl_data = self._load_smiles_file(self.chemspace_path / "ertl_data.smi")
        logger.info(f"Loaded {len(self.ertl_data)} molecules from ertl_data.smi")
        
        # Load chemspace data (medium priority)
        self.chemspace_data = self._load_chemspace_data()
        logger.info(f"Loaded {len(self.chemspace_data)} molecules from chemspace collections")
        
        # Load cluster centroids (packed fingerprints) (lowest priority)
        centroids_file = self.cluster_path / "cluster-centroids-packed.pkl"
        if not centroids_file.exists():
            raise FileNotFoundError(f"Cluster centroids file not found: {centroids_file}")
        
        with open(centroids_file, "rb") as f:
            self.centroids = pickle.load(f)
        
        logger.info(f"Loaded {len(self.centroids)} cluster centroids")
        
        # Load ring SMILES for mapping cluster IDs to actual structures
        self.ring_smiles = self._load_ring_smiles()
        
        # Load cluster assignments to map centroids to ring IDs
        clusters_file = self.cluster_path / "clusters.pkl"
        if clusters_file.exists():
            with open(clusters_file, "rb") as f:
                self.clusters = pickle.load(f)
            logger.info(f"Loaded cluster assignments for {len(self.clusters)} clusters")
        else:
            self.clusters = None
            logger.warning("Cluster assignments file not found - will only return centroid matches")
        
        # Load fingerprints if available for exact lookups
        fps_file = self.cluster_path / "input-fps" / "rings.npy"
        if fps_file.exists():
            self.fingerprints = np.load(fps_file)
            logger.info(f"Loaded {len(self.fingerprints)} fingerprints")
        else:
            self.fingerprints = None
    
    def _load_smiles_file(self, filepath: Path) -> List[str]:
        """
        Load SMILES from a .smi file.
        
        Args:
            filepath: Path to .smi file
            
        Returns:
            List of SMILES strings
        """
        if not filepath.exists():
            logger.warning(f"SMILES file not found: {filepath}")
            return []
        
        smiles_list = []
        with open(filepath, "r") as f:
            for line in f:
                line = line.strip()
                # Skip comments and empty lines
                if line.startswith("#") or not line:
                    continue
                # SMILES is the first column (tab-separated)
                parts = line.split("\t")
                if parts and parts[0]:
                    smiles_list.append(parts[0])
        
        return smiles_list
    
    def _load_chemspace_data(self) -> List[str]:
        """Load all chemspace .smi files (excluding ertl_data)."""
        all_smiles = []
        
        smi_files = [
            "Chemspace_Bioactive_Rings_Focused_Collection.smi",
            "Chemspace_Carboxylic_Acid_Bioisosteres.smi",
            "Chemspace_meta-Substituted_Benzene_Bioisosteres.smi",
            "Chemspace_ortho-Substituted_Benzene_Bioisosteres.smi",
            "Chemspace_para-Substituted_Benzene_Bioisosteres.smi",
            "Chemspace_tert-Butyl_Group_Bioisosteres.smi",
        ]
        
        for smi_file in smi_files:
            filepath = self.chemspace_path / smi_file
            smiles = self._load_smiles_file(filepath)
            all_smiles.extend(smiles)
            logger.info(f"  Loaded {len(smiles)} from {smi_file}")
        
        return all_smiles
    
    def _load_ring_smiles(self) -> List[str]:
        """Load ring SMILES from the database file."""
        rings_file = self.db_path / "rings.smi"
        if not rings_file.exists():
            raise FileNotFoundError(f"Ring SMILES file not found: {rings_file}")
        
        smiles_list = []
        with open(rings_file, "r") as f:
            for line in f:
                line = line.strip()
                # Skip comments and empty lines
                if line.startswith("#") or not line:
                    continue
                # SMILES is the first column (tab-separated)
                parts = line.split("\t")
                if parts:
                    smiles_list.append(parts[0])
        
        logger.info(f"Loaded {len(smiles_list)} ring SMILES")
        return smiles_list
    
    def _calculate_descriptors(self, smiles: str) -> dict:
        """
        Calculate physicochemical descriptors for bio-isostere comparison.
        
        Args:
            smiles: SMILES string of the ring
            
        Returns:
            Dictionary of descriptors including pharmacophore and ring topology
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {}
        
        try:
            # Get ring information
            ring_info = mol.GetRingInfo()
            ring_sizes = [len(ring) for ring in ring_info.AtomRings()]
            
            # Count fused rings (rings sharing atoms)
            num_fused = 0
            atom_rings = ring_info.AtomRings()
            for i, ring1 in enumerate(atom_rings):
                for ring2 in atom_rings[i+1:]:
                    if len(set(ring1) & set(ring2)) > 0:
                        num_fused += 1
            
            descriptors = {
                # Lipophilicity
                'logp': Crippen.MolLogP(mol),
                
                # Polarity / H-bonding
                'tpsa': Descriptors.TPSA(mol),
                'h_donors': Lipinski.NumHDonors(mol),
                'h_acceptors': Lipinski.NumHAcceptors(mol),
                
                # Aromaticity
                'aromatic_rings': Descriptors.NumAromaticRings(mol),
                'aromatic_atoms': len([a for a in mol.GetAtoms() if a.GetIsAromatic()]),
                
                # Heteroatoms
                'num_heteroatoms': Lipinski.NumHeteroatoms(mol),
                'num_n': len([a for a in mol.GetAtoms() if a.GetAtomicNum() == 7]),
                'num_o': len([a for a in mol.GetAtoms() if a.GetAtomicNum() == 8]),
                'num_s': len([a for a in mol.GetAtoms() if a.GetAtomicNum() == 16]),
                
                # Size / complexity
                'mw': Descriptors.MolWt(mol),
                'num_rings': Descriptors.RingCount(mol),
                'num_atoms': mol.GetNumAtoms(),
                'num_heavy_atoms': mol.GetNumHeavyAtoms(),
                
                # Ring topology
                'ring_sizes': sorted(ring_sizes) if ring_sizes else [],
                'num_fused_rings': num_fused,
                'largest_ring': max(ring_sizes) if ring_sizes else 0,
                'smallest_ring': min(ring_sizes) if ring_sizes else 0,
            }
            return descriptors
        except Exception as e:
            logger.error(f"Error calculating descriptors for {smiles}: {e}")
            return {}
    
    def _calculate_2d_pharmacophore_similarity(self, query_smiles: str, match_smiles: str) -> float:
        """
        Calculate 2D pharmacophore similarity using Gobbi pharmacophore fingerprints.
        
        Args:
            query_smiles: Query molecule SMILES
            match_smiles: Match molecule SMILES
            
        Returns:
            Pharmacophore Tanimoto similarity (0-1)
        """
        try:
            query_mol = Chem.MolFromSmiles(query_smiles)
            match_mol = Chem.MolFromSmiles(match_smiles)
            
            if query_mol is None or match_mol is None:
                return 0.0
            
            # Generate Gobbi 2D pharmacophore fingerprints
            factory = Gobbi_Pharm2D.factory
            query_fp = Generate.Gen2DFingerprint(query_mol, factory)
            match_fp = Generate.Gen2DFingerprint(match_mol, factory)
            
            # Calculate Tanimoto similarity
            similarity = Chem.DataStructs.TanimotoSimilarity(query_fp, match_fp)
            return float(similarity)
            
        except Exception as e:
            logger.error(f"Error calculating pharmacophore similarity: {e}")
            return 0.0
    
    def _calculate_ring_topology_similarity(self, query_desc: dict, match_desc: dict) -> float:
        """
        Calculate ring topology similarity based on ring sizes and fusion patterns.
        
        Strict matching to ensure single rings match single rings, not fused systems.
        
        Args:
            query_desc: Query molecule descriptors
            match_desc: Match molecule descriptors
            
        Returns:
            Ring topology similarity score (0-1)
        """
        try:
            query_rings = query_desc.get('ring_sizes', [])
            match_rings = match_desc.get('ring_sizes', [])
            
            if not query_rings or not match_rings:
                return 0.0
            
            # Exact match of ring sizes gets score of 1.0
            if query_rings == match_rings:
                return 1.0
            
            # STRICT: If number of rings differs, apply heavy penalty
            # Single rings should NOT match fused ring systems
            query_num_rings = len(query_rings)
            match_num_rings = len(match_rings)
            
            if query_num_rings != match_num_rings:
                # Different number of rings = maximum 0.3 similarity
                # This prevents single rings from matching multi-ring systems
                ring_count_penalty = abs(query_num_rings - match_num_rings)
                if ring_count_penalty >= 2:
                    return 0.0  # 2+ ring difference = incompatible
                else:
                    # 1 ring difference = max 0.3 score
                    max_score = 0.3
            else:
                # Same number of rings = can score up to 1.0
                max_score = 1.0
            
            # Ring size similarity (Jaccard)
            ring_intersection = len(set(query_rings) & set(match_rings))
            ring_union = len(set(query_rings) | set(match_rings))
            ring_size_sim = ring_intersection / ring_union if ring_union > 0 else 0.0
            
            # Fusion pattern similarity (only matters for multi-ring systems)
            query_fused = query_desc.get('num_fused_rings', 0)
            match_fused = match_desc.get('num_fused_rings', 0)
            
            if query_num_rings == 1 and match_num_rings == 1:
                # Both single rings - fusion doesn't matter
                fusion_sim = 1.0
            elif query_fused == 0 and match_fused == 0:
                # Both have no fusion (separate rings)
                fusion_sim = 1.0
            elif query_fused == 0 or match_fused == 0:
                # One fused, one not = incompatible
                fusion_sim = 0.0
            else:
                # Both fused - check similarity
                fusion_sim = 1.0 - min(abs(query_fused - match_fused) / max(query_fused, match_fused), 1.0)
            
            # Weighted combination
            topology_sim = (
                0.5 * ring_size_sim +      # 50% ring sizes
                0.5 * fusion_sim            # 50% fusion pattern
            )
            
            # Apply max score cap
            topology_sim = min(topology_sim, max_score)
            
            return float(topology_sim)
            
        except Exception as e:
            logger.error(f"Error calculating ring topology similarity: {e}")
            return 0.0
    
    def _check_bio_isostere_filters(
        self,
        query_desc: dict,
        match_desc: dict,
        logp_tolerance: float = 1.5,
        tpsa_tolerance: float = 25.0,
        h_donor_tolerance: int = 1,
        h_acceptor_tolerance: int = 2
    ) -> Tuple[bool, dict]:
        """
        Check if a match passes bio-isostere filters.
        
        Args:
            query_desc: Query molecule descriptors
            match_desc: Match molecule descriptors
            logp_tolerance: Maximum logP difference
            tpsa_tolerance: Maximum TPSA difference (Ų)
            h_donor_tolerance: Maximum H-donor count difference
            h_acceptor_tolerance: Maximum H-acceptor count difference
            
        Returns:
            Tuple of (passes_filters, differences_dict)
        """
        if not query_desc or not match_desc:
            return False, {}
        
        # Calculate differences
        diffs = {
            'delta_logp': abs(query_desc.get('logp', 0) - match_desc.get('logp', 0)),
            'delta_tpsa': abs(query_desc.get('tpsa', 0) - match_desc.get('tpsa', 0)),
            'delta_h_donors': abs(query_desc.get('h_donors', 0) - match_desc.get('h_donors', 0)),
            'delta_h_acceptors': abs(query_desc.get('h_acceptors', 0) - match_desc.get('h_acceptors', 0)),
            'delta_aromatic_rings': abs(query_desc.get('aromatic_rings', 0) - match_desc.get('aromatic_rings', 0)),
            'delta_heteroatoms': abs(query_desc.get('num_heteroatoms', 0) - match_desc.get('num_heteroatoms', 0)),
        }
        
        # Check filters
        passes = (
            diffs['delta_logp'] <= logp_tolerance and
            diffs['delta_tpsa'] <= tpsa_tolerance and
            diffs['delta_h_donors'] <= h_donor_tolerance and
            diffs['delta_h_acceptors'] <= h_acceptor_tolerance
        )
        
        return passes, diffs
    
    def _smiles_to_packed_fingerprint(self, smiles: str, fp_size: int = 2048) -> np.ndarray:
        """
        Convert SMILES to packed Morgan fingerprint.
        
        Args:
            smiles: SMILES string of the ring
            fp_size: Fingerprint size (default 2048)
            
        Returns:
            Packed fingerprint as numpy array
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {smiles}")
        
        # Generate Morgan fingerprint (radius 2, 2048 bits)
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=fp_size)
        
        # Convert to numpy array
        fp_array = np.zeros((fp_size,), dtype=np.uint8)
        for i in range(fp_size):
            fp_array[i] = fp[i]
        
        # Pack the fingerprint (8 bits per byte)
        n_bytes = (fp_size + 7) // 8
        packed_fp = np.zeros(n_bytes, dtype=np.uint8)
        
        for i in range(fp_size):
            if fp_array[i]:
                byte_idx = i // 8
                bit_idx = 7 - (i % 8)
                packed_fp[byte_idx] |= (1 << bit_idx)
        
        return packed_fp
    
    def _calculate_similarities(self, query_fp: np.ndarray) -> np.ndarray:
        """
        Calculate Tanimoto similarities between query fingerprint and all centroids.
        
        Args:
            query_fp: Packed query fingerprint
            
        Returns:
            Array of similarity scores
        """
        # Calculate cardinality of query fingerprint (number of set bits)
        query_card = np.unpackbits(query_fp).sum()
        
        # Calculate cardinalities of all centroids
        centroid_cards = np.array([np.unpackbits(c).sum() for c in self.centroids])
        
        # Calculate similarities using bblean's optimized function
        similarities = np.zeros(len(self.centroids), dtype=np.float32)
        
        for i, centroid in enumerate(self.centroids):
            # Calculate intersection (AND operation)
            intersection = np.bitwise_and(query_fp, centroid)
            intersection_card = np.unpackbits(intersection).sum()
            
            # Tanimoto similarity = intersection / union
            union_card = query_card + centroid_cards[i] - intersection_card
            if union_card > 0:
                similarities[i] = intersection_card / union_card
            else:
                similarities[i] = 0.0
        
        return similarities
    
    
    def scan(
        self,
        query_smiles: str,
        top_k: int = 20,
        min_similarity: float = 0.3,
        use_bio_filters: bool = True,
        logp_tolerance: float = 1.5,
        tpsa_tolerance: float = 25.0,
        h_donor_tolerance: int = 1,
        h_acceptor_tolerance: int = 2
    ) -> List[dict]:
        """
        Hierarchical scan for bio-isosteric replacements.
        
        Scans in priority order:
        1. ertl_data.smi (highest priority)
        2. chemspace collections (medium priority)
        3. cluster centroids (lowest priority)
        
        Only replaces candidates if later source has BETTER bio-isostere score.
        
        Args:
            query_smiles: SMILES string of the query ring
            top_k: Number of top matches to return (default 20)
            min_similarity: Minimum Tanimoto similarity threshold (default 0.3)
            use_bio_filters: Apply bio-isostere physicochemical filters (default True)
            logp_tolerance: Maximum logP difference (default 1.5)
            tpsa_tolerance: Maximum TPSA difference in Ų (default 25.0)
            h_donor_tolerance: Maximum H-donor count difference (default 1)
            h_acceptor_tolerance: Maximum H-acceptor count difference (default 2)
            
        Returns:
            List of dictionaries containing match information sorted by bio-isostere score
        """
        # Calculate descriptors for query if bio-filters are enabled
        query_desc = {}
        if use_bio_filters:
            query_desc = self._calculate_descriptors(query_smiles)
            if not query_desc:
                logger.warning("Could not calculate descriptors for query - disabling bio-filters")
                use_bio_filters = False
        
        # Dictionary to store best match per canonical SMILES
        # Key: canonical SMILES, Value: match dict with bio_isostere_score
        best_matches = {}
        
        # Exclude query itself
        query_canonical = None
        try:
            query_mol = Chem.MolFromSmiles(query_smiles)
            if query_mol is not None:
                query_canonical = Chem.MolToSmiles(query_mol)
        except Exception:
            pass
        
        # Step 1: Scan ertl_data (highest priority)
        logger.info(f"Scanning ertl_data ({len(self.ertl_data)} molecules)...")
        self._scan_smiles_list(
            query_smiles, query_desc, self.ertl_data, 
            min_similarity, use_bio_filters,
            logp_tolerance, tpsa_tolerance, h_donor_tolerance, h_acceptor_tolerance,
            best_matches, query_canonical, source="ertl"
        )
        
        # Step 2: Scan chemspace (medium priority)
        logger.info(f"Scanning chemspace ({len(self.chemspace_data)} molecules)...")
        self._scan_smiles_list(
            query_smiles, query_desc, self.chemspace_data,
            min_similarity, use_bio_filters,
            logp_tolerance, tpsa_tolerance, h_donor_tolerance, h_acceptor_tolerance,
            best_matches, query_canonical, source="chemspace"
        )
        
        # Step 3: Scan clusters (lowest priority)
        logger.info(f"Scanning clusters ({len(self.centroids)} centroids)...")
        self._scan_clusters(
            query_smiles, query_desc,
            min_similarity, use_bio_filters,
            logp_tolerance, tpsa_tolerance, h_donor_tolerance, h_acceptor_tolerance,
            best_matches, query_canonical
        )
        
        # Convert to list and sort by bio-isostere score
        results = list(best_matches.values())
        results.sort(key=lambda x: x['bio_isostere_score'], reverse=True)
        
        logger.info(f"Found {len(results)} unique bio-isosteres")
        
        # Return top_k results
        return results[:top_k]
    
    def _scan_smiles_list(
        self,
        query_smiles: str,
        query_desc: dict,
        smiles_list: List[str],
        min_similarity: float,
        use_bio_filters: bool,
        logp_tolerance: float,
        tpsa_tolerance: float,
        h_donor_tolerance: int,
        h_acceptor_tolerance: int,
        best_matches: dict,
        query_canonical: Optional[str],
        source: str
    ):
        """Scan a list of SMILES and update best_matches."""
        # Use new MorganGenerator API to avoid deprecation warnings
        morgan_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
        query_mol = Chem.MolFromSmiles(query_smiles)
        query_fp = morgan_gen.GetFingerprint(query_mol)
        
        for smiles in smiles_list:
            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    continue
                
                canonical = Chem.MolToSmiles(mol)
                
                # Skip query itself
                if query_canonical and canonical == query_canonical:
                    continue
                
                # Calculate similarity using new API
                match_fp = morgan_gen.GetFingerprint(mol)
                similarity = Chem.DataStructs.TanimotoSimilarity(query_fp, match_fp)
                
                if similarity < min_similarity:
                    continue
                
                # Calculate bio-isostere score
                match_desc = self._calculate_descriptors(smiles)
                if not match_desc:
                    continue
                
                # Calculate property differences (for information only, no filtering)
                if use_bio_filters:
                    _, delta_props = self._check_bio_isostere_filters(
                        query_desc, match_desc,
                        logp_tolerance, tpsa_tolerance,
                        h_donor_tolerance, h_acceptor_tolerance
                    )
                else:
                    delta_props = {}
                
                # Calculate component scores
                pharmacophore_sim = self._calculate_2d_pharmacophore_similarity(
                    query_smiles, smiles
                )
                topology_sim = self._calculate_ring_topology_similarity(
                    query_desc, match_desc
                ) if use_bio_filters else 0.0
                
                # Adaptive weighting based on query molecule size
                query_heavy_atoms = query_desc.get('num_heavy_atoms', 20)
                if query_heavy_atoms < 10:
                    # Small molecules: trust pharmacophore/topology more
                    # Morgan FP similarity is unreliable for small molecules
                    w_morgan, w_pharm, w_topo = 0.2, 0.5, 0.3
                else:
                    # Large molecules: Morgan FP more reliable
                    w_morgan, w_pharm, w_topo = 0.4, 0.4, 0.2
                
                bio_isostere_score = (
                    w_morgan * similarity +
                    w_pharm * pharmacophore_sim +
                    w_topo * topology_sim
                )
                
                # Only add/update if better than existing
                if canonical not in best_matches or bio_isostere_score > best_matches[canonical]['bio_isostere_score']:
                    best_matches[canonical] = {
                        "source": source,
                        "similarity": float(similarity),
                        "centroid_smiles": smiles,
                        "bio_isostere_score": float(bio_isostere_score),
                        "pharmacophore_similarity": float(pharmacophore_sim),
                        "topology_similarity": float(topology_sim),
                        "descriptors": match_desc,
                        "delta_properties": delta_props
                    }
                    
            except Exception as e:
                logger.debug(f"Error processing {smiles}: {e}")
                continue
    
    def _scan_clusters(
        self,
        query_smiles: str,
        query_desc: dict,
        min_similarity: float,
        use_bio_filters: bool,
        logp_tolerance: float,
        tpsa_tolerance: float,
        h_donor_tolerance: int,
        h_acceptor_tolerance: int,
        best_matches: dict,
        query_canonical: Optional[str]
    ):
        """Scan cluster centroids and update best_matches."""
        # Generate fingerprint for query
        query_fp = self._smiles_to_packed_fingerprint(query_smiles)
        
        # Calculate similarities
        similarities = self._calculate_similarities(query_fp)
        
        # Filter by minimum similarity
        valid_indices = np.where(similarities >= min_similarity)[0]
        
        if len(valid_indices) == 0:
            return
        
        # Sort by similarity (descending) - get many candidates
        sorted_indices = valid_indices[np.argsort(-similarities[valid_indices])]
        
        for cluster_id in sorted_indices:
            similarity = float(similarities[cluster_id])
            
            # Get cluster members if available
            num_members = 0
            example_smiles = []
            
            if self.clusters is not None and cluster_id < len(self.clusters):
                member_indices = self.clusters[cluster_id]
                num_members = len(member_indices)
                
                # Get up to 5 example SMILES from this cluster
                example_indices = member_indices[:5] if len(member_indices) > 5 else member_indices
                example_smiles = [
                    self.ring_smiles[idx] 
                    for idx in example_indices 
                    if idx < len(self.ring_smiles)
                ]
            
            # The centroid itself represents the cluster
            # We'll use the first member as the representative if available
            centroid_smiles = example_smiles[0] if example_smiles else "N/A"
            
            # Skip if N/A
            if centroid_smiles == "N/A":
                continue
                
            try:
                mol = Chem.MolFromSmiles(centroid_smiles)
                if mol is None:
                    continue
                
                canonical = Chem.MolToSmiles(mol)
                
                # Skip query itself
                if query_canonical and canonical == query_canonical:
                    continue
                
            except Exception:
                continue
            
            # Calculate descriptors and apply bio-isostere filters
            match_desc = {}
            delta_props = {}
            bio_isostere_score = similarity  # Default to structural similarity
            passes_filters = True
            pharmacophore_sim = 0.0
            topology_sim = 0.0
            
            if use_bio_filters and centroid_smiles != "N/A":
                match_desc = self._calculate_descriptors(centroid_smiles)
                if match_desc:
                    # Calculate property differences (for information only)
                    _, delta_props = self._check_bio_isostere_filters(
                        query_desc,
                        match_desc,
                        logp_tolerance=logp_tolerance,
                        tpsa_tolerance=tpsa_tolerance,
                        h_donor_tolerance=h_donor_tolerance,
                        h_acceptor_tolerance=h_acceptor_tolerance
                    )
                    
                    # Calculate 2D pharmacophore similarity
                    pharmacophore_sim = self._calculate_2d_pharmacophore_similarity(
                        query_smiles,
                        centroid_smiles
                    )
                    
                    # Calculate ring topology similarity
                    topology_sim = self._calculate_ring_topology_similarity(
                        query_desc,
                        match_desc
                    )
                    
                    # Adaptive weighting based on query molecule size
                    query_heavy_atoms = query_desc.get('num_heavy_atoms', 20)
                    if query_heavy_atoms < 10:
                        # Small molecules: trust pharmacophore/topology more
                        # Morgan FP similarity is unreliable for small molecules
                        w_morgan, w_pharm, w_topo = 0.2, 0.5, 0.3
                    else:
                        # Large molecules: Morgan FP more reliable
                        w_morgan, w_pharm, w_topo = 0.4, 0.4, 0.2
                    
                    # Calculate composite bio-isostere score with adaptive weights
                    bio_isostere_score = (
                        w_morgan * similarity +
                        w_pharm * pharmacophore_sim +
                        w_topo * topology_sim
                    )
            
            # Only add/update if better than existing
            if canonical not in best_matches or bio_isostere_score > best_matches[canonical]['bio_isostere_score']:
                best_matches[canonical] = {
                    "source": "clusters",
                    "cluster_id": int(cluster_id),
                    "similarity": float(similarity),
                    "centroid_smiles": centroid_smiles,
                    "num_members": num_members,
                    "example_smiles": example_smiles,
                    "bio_isostere_score": float(bio_isostere_score),
                    "pharmacophore_similarity": float(pharmacophore_sim),
                    "topology_similarity": float(topology_sim),
                    "descriptors": match_desc,
                    "delta_properties": delta_props
                }
    
    def get_cluster_info(self, cluster_id: int) -> Optional[dict]:
        """
        Get detailed information about a specific cluster.
        
        Args:
            cluster_id: ID of the cluster
            
        Returns:
            Dictionary with cluster information or None if not found
        """
        if cluster_id >= len(self.centroids):
            return None
        
        info = {
            "cluster_id": cluster_id,
            "num_members": 0,
            "member_smiles": []
        }
        
        if self.clusters is not None and cluster_id < len(self.clusters):
            member_indices = self.clusters[cluster_id]
            info["num_members"] = len(member_indices)
            info["member_smiles"] = [
                self.ring_smiles[idx] 
                for idx in member_indices 
                if idx < len(self.ring_smiles)
            ]
        
        return info


# Singleton instance
_scanner_instance = None


def get_scanner() -> BioisostereScanner:
    """Get or create the singleton BioisostereScanner instance."""
    global _scanner_instance
    if _scanner_instance is None:
        _scanner_instance = BioisostereScanner()
    return _scanner_instance


if __name__ == "__main__":
    scanner = get_scanner()
    print(scanner.scan("c1ccccc1"))