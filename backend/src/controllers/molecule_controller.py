"""
Controller for molecule analysis operations.

Handles the business logic for analyzing molecules using SMILES strings
and interfacing with the SmartChemist tool.
"""
from typing import List, Tuple, Optional, Dict, Any
import math
import os
import re
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, rdMolTransforms, rdMolAlign, rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Geometry import Point3D
import numpy as np
import copy
from openai import OpenAI

from src.services.mcp.smart_chemist.tools.smart_chemist import SmartChemist
from src.api.models.molecule import MoleculeResponse, Match, TrivialName, HighlightResponse


class MoleculeController:
    """Controller for molecule analysis operations."""
    
    def __init__(self):
        """Initialize the MoleculeController with a SmartChemist instance."""
        self.smart_chemist = SmartChemist()
    
    def analyze_molecule(self, smiles: str) -> MoleculeResponse:
        """
        Analyze a molecule from its SMILES string.
        
        Args:
            smiles: SMILES string representation of the molecule
            
        Returns:
            MoleculeResponse containing the annotated molecule information
            
        Raises:
            ValueError: If the SMILES string is invalid
        """
        # Convert SMILES to RDKit molecule object
        mol = Chem.MolFromSmiles(smiles)
        
        if mol is None:
            raise ValueError(f"Invalid SMILES string: {smiles}")
        
        # Get annotation JSON from SmartChemist
        annotation_data = self.smart_chemist.mol_to_annotation_json(mol)
        
        # Convert matches to the proper format
        matches = []
        for match_data in annotation_data["matches"]:
            trivial_name = TrivialName(
                name=match_data["trivial_name"]["name"],
                smarts=match_data["trivial_name"]["smarts"],
                group=match_data["trivial_name"]["group"],
                bonds=match_data["trivial_name"]["bonds"],
                hierarchy=match_data["trivial_name"]["hierarchy"],
                index=match_data["trivial_name"]["index"]
            )
            match = Match(
                atom_indices=list(match_data["atom_indices"]),
                svg=match_data.get("svg"),
                trivial_name=trivial_name
            )
            matches.append(match)
        
        # Create and return the response
        return MoleculeResponse(
            name=annotation_data["name"],
            svg=annotation_data["svg"],
            matches=matches,
            smiles=annotation_data["smiles"]
        )
    
    def highlight_molecule(self, smiles: str, atom_indices: List[int]) -> HighlightResponse:
        """
        Generate an SVG of a molecule with specific atoms highlighted.
        
        Args:
            smiles: SMILES string representation of the molecule
            atom_indices: List of atom indices to highlight
            
        Returns:
            HighlightResponse containing the highlighted SVG
            
        Raises:
            ValueError: If the SMILES string is invalid
        """
        mol = Chem.MolFromSmiles(smiles)
        
        if mol is None:
            raise ValueError(f"Invalid SMILES string: {smiles}")
        
        # Generate 2D coordinates
        AllChem.Compute2DCoords(mol)
        
        # Create the drawer with the same size as the source molecule
        drawer = rdMolDraw2D.MolDraw2DSVG(400, 400)
        
        # Set transparent background
        drawer.drawOptions().clearBackground = False
        
        # Set highlight color (light blue/cyan)
        highlight_color = (0.7, 0.9, 1.0)  # Light blue
        
        # Create highlight atom and bond dictionaries
        highlight_atoms = {idx: highlight_color for idx in atom_indices}
        
        # Find bonds between highlighted atoms
        highlight_bonds = {}
        for bond in mol.GetBonds():
            begin_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            if begin_idx in atom_indices and end_idx in atom_indices:
                highlight_bonds[bond.GetIdx()] = highlight_color
        
        # Draw molecule with highlights
        drawer.DrawMolecule(
            mol,
            highlightAtoms=list(atom_indices),
            highlightAtomColors=highlight_atoms,
            highlightBonds=list(highlight_bonds.keys()),
            highlightBondColors=highlight_bonds
        )
        drawer.FinishDrawing()
        
        svg = drawer.GetDrawingText()
        
        return HighlightResponse(
            svg=svg,
            smiles=Chem.MolToSmiles(mol)
        )
    
    def extract_substructure_smiles(self, smiles: str, atom_indices: List[int]) -> str:
        """
        Extract a substructure SMILES from a molecule given atom indices.
        
        Args:
            smiles: SMILES string representation of the parent molecule
            atom_indices: List of atom indices that define the substructure
            
        Returns:
            SMILES string of the extracted substructure
            
        Raises:
            ValueError: If the SMILES string is invalid or atom indices are invalid
        """
        mol = Chem.MolFromSmiles(smiles)
        
        if mol is None:
            raise ValueError(f"Invalid SMILES string: {smiles}")
        
        if not atom_indices:
            raise ValueError("No atom indices provided")
        
        # Validate atom indices
        num_atoms = mol.GetNumAtoms()
        if any(idx < 0 or idx >= num_atoms for idx in atom_indices):
            raise ValueError(f"Invalid atom indices. Molecule has {num_atoms} atoms")
        
        # Create a new molecule with only the selected atoms
        # First, get the bonds between the selected atoms
        atom_set = set(atom_indices)
        bonds_to_include = []
        
        for bond in mol.GetBonds():
            begin_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            if begin_idx in atom_set and end_idx in atom_set:
                bonds_to_include.append(bond.GetIdx())
        
        # Use PathToSubmol to extract the substructure
        # This preserves the connectivity and chemistry
        submol = Chem.PathToSubmol(mol, bonds_to_include, atomMap=atom_indices)
        
        if submol is None:
            # Fallback: try to create a fragment
            submol = Chem.RWMol()
            atom_map = {}
            
            # Add atoms
            for idx in atom_indices:
                atom = mol.GetAtomWithIdx(idx)
                new_idx = submol.AddAtom(Chem.Atom(atom.GetSymbol()))
                atom_map[idx] = new_idx
            
            # Add bonds between the atoms
            for bond in mol.GetBonds():
                begin_idx = bond.GetBeginAtomIdx()
                end_idx = bond.GetEndAtomIdx()
                if begin_idx in atom_set and end_idx in atom_set:
                    submol.AddBond(
                        atom_map[begin_idx],
                        atom_map[end_idx],
                        bond.GetBondType()
                    )
            
            submol = submol.GetMol()
        
        if submol is None:
            raise ValueError("Failed to extract substructure")
        
        # Sanitize and return SMILES
        try:
            Chem.SanitizeMol(submol)
            return Chem.MolToSmiles(submol)
        except Exception as e:
            # If sanitization fails, return without sanitization
            return Chem.MolToSmiles(submol, sanitize=False)
    
    def replace_substructure_with_alignment(
        self, 
        parent_smiles: str, 
        source_pattern_smiles: str,
        replacement_smiles: str,
        atom_indices: List[int]
    ) -> List[str]:
        """
        Replace a substructure in a molecule with a replacement pattern.
        Fully connection-safe, multi-attempt ring replacement with RMSD scoring.
        Manually builds the molecule to ensure proper connectivity.
        
        Args:
            parent_smiles: SMILES of the parent molecule
            source_pattern_smiles: SMILES of the pattern to replace
            replacement_smiles: SMILES of the replacement pattern (bio-isostere)
            atom_indices: Indices of atoms in the source pattern within the parent molecule
            
        Returns:
            List of SMILES strings for successfully generated molecules (up to 3)
            
        Raises:
            ValueError: If inputs are invalid
        """
        print(f"\n=== Replace Substructure Debug ===")
        print(f"Parent: {parent_smiles}")
        print(f"Source pattern: {source_pattern_smiles}")
        print(f"Replacement: {replacement_smiles}")
        print(f"Atom indices: {atom_indices}")
        
        parent_mol = Chem.MolFromSmiles(parent_smiles)
        old_ring = Chem.MolFromSmiles(source_pattern_smiles)
        new_ring = Chem.MolFromSmiles(replacement_smiles)
        
        if not parent_mol or not old_ring or not new_ring:
            raise ValueError("Invalid SMILES string(s)")
        
        # Generate 2D coordinates
        rdDepictor.Compute2DCoords(parent_mol)
        rdDepictor.Compute2DCoords(old_ring)
        rdDepictor.Compute2DCoords(new_ring)
        
        try:
            # Find matching occurrences of old ring in parent
            matches = parent_mol.GetSubstructMatches(old_ring)
            
            if not matches:
                print("Old ring not found in parent molecule")
                return [parent_smiles]
            
            # Use the atom_indices if they match, otherwise use first match
            old_match = None
            for match in matches:
                if set(match) == set(atom_indices):
                    old_match = match
                    break
            
            if old_match is None:
                old_match = matches[0]
            
            old_match_set = set(old_match)
            
            # Find exit vectors (bonds leaving the ring system)
            exit_vectors = []  # list of (old_ring_atom, external_atom, BondType)
            for old_idx in old_match:
                atom = parent_mol.GetAtomWithIdx(old_idx)
                for nbr in atom.GetNeighbors():
                    nbr_idx = nbr.GetIdx()
                    if nbr_idx not in old_match_set:
                        bond = parent_mol.GetBondBetweenAtoms(old_idx, nbr_idx)
                        exit_vectors.append((old_idx, nbr_idx, bond.GetBondType()))
            
            print(f"Exit vectors found: {len(exit_vectors)}")
            for ev in exit_vectors:
                print(f"  {ev[0]} (ring) -> {ev[1]} (external), bond: {ev[2]}")
            
            if not exit_vectors:
                print("Selected old ring has no exit vectors; cannot replace.")
                return [parent_smiles]
            
            # Collect all successful results
            successful_results = []
            
            # Try 3 attempts with different anchor strategies
            attempts = 3
            for trial in range(attempts):
                try:
                    print(f"\n--- Attempt {trial} ---")
                    # Generate anchor mapping
                    n_old = old_ring.GetNumAtoms()
                    n_new = new_ring.GetNumAtoms()
                    k = min(2, min(n_old, n_new))
                    
                    np.random.seed(42 + trial * 100)
                    chosen_old = np.random.choice(range(n_old), k, replace=False)
                    chosen_new = np.random.choice(range(n_new), k, replace=False)
                    atom_map = list(zip(chosen_new, chosen_old))
                    print(f"Atom map: {atom_map}")
                    
                    # Align new ring to old ring
                    new_aligned = Chem.Mol(new_ring)
                    try:
                        rdDepictor.GenerateDepictionMatching2DStructure(
                            new_aligned, old_ring, atomMap=atom_map
                        )
                        print("✓ Alignment succeeded")
                    except Exception as e:
                        print(f"✗ Alignment failed: {str(e)}")
                        continue
                    
                    # Build new molecule manually (safe)
                    em = Chem.EditableMol(Chem.Mol())
                    
                    # Map: new_idx → old_idx
                    new_to_old = dict(atom_map)
                    
                    # Add atoms from parent (except old ring)
                    parent_to_newidx = {}
                    for i, atom in enumerate(parent_mol.GetAtoms()):
                        if i not in old_match_set:
                            new_i = em.AddAtom(copy.copy(atom))
                            parent_to_newidx[i] = new_i
                    
                    # Add all atoms of new ring
                    newring_to_newidx = {}
                    for j, atom in enumerate(new_aligned.GetAtoms()):
                        new_idx = em.AddAtom(copy.copy(atom))
                        newring_to_newidx[j] = new_idx
                    
                    # Add bonds for parent (except bonds inside old ring)
                    for bond in parent_mol.GetBonds():
                        a1, a2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
                        if a1 in old_match_set or a2 in old_match_set:
                            continue
                        em.AddBond(
                            parent_to_newidx[a1],
                            parent_to_newidx[a2],
                            bond.GetBondType()
                        )
                    
                    # Add bonds inside the new ring
                    for bond in new_aligned.GetBonds():
                        a1 = newring_to_newidx[bond.GetBeginAtomIdx()]
                        a2 = newring_to_newidx[bond.GetEndAtomIdx()]
                        em.AddBond(a1, a2, bond.GetBondType())
                    
                    # Reconnect exit vectors
                    print(f"  Reconnecting {len(exit_vectors)} exit vectors...")
                    valid = True
                    for (old_anchor, ext_atom, bondtype) in exit_vectors:
                        # Find which old ring atom index this is
                        old_ring_idx = old_match.index(old_anchor) if old_anchor in old_match else None
                        if old_ring_idx is None:
                            print(f"    ✗ old_anchor {old_anchor} not found in old_match {old_match}")
                            valid = False
                            break
                        
                        # Find new ring atom that corresponds
                        new_candidates = [new_a for (new_a, old_a) in atom_map if old_a == old_ring_idx]
                        
                        if not new_candidates:
                            print(f"    ✗ No new ring atom maps to old_ring_idx {old_ring_idx}, atom_map={atom_map}")
                            valid = False
                            break
                        
                        new_anchor = newring_to_newidx[new_candidates[0]]
                        ext_new = parent_to_newidx[ext_atom]
                        
                        print(f"    Connecting: new_ring[{new_candidates[0]}]={new_anchor} -> parent[{ext_atom}]={ext_new}, bond={bondtype}")
                        em.AddBond(new_anchor, ext_new, bondtype)
                    
                    if not valid:
                        print(f"  ✗ Exit vector reconnection validation failed")
                        continue
                    
                    print("  ✓ All exit vectors reconnected")
                    
                    # Final molecule
                    candidate = em.GetMol()
                    print(f"  Built molecule with {candidate.GetNumAtoms()} atoms")
                    
                    # Sanitize
                    try:
                        Chem.SanitizeMol(candidate)
                        print("  ✓ Sanitization succeeded")
                    except Exception as e:
                        print(f"  ✗ Sanitization failed: {str(e)}")
                        continue
                    
                    if not self._is_connected(candidate):
                        print(f"  ✗ Molecule not connected ({len(Chem.GetMolFrags(candidate))} fragments)")
                        continue
                    
                    print("  ✓ Molecule is connected")
                    
                    # Generate 2D coords and get SMILES
                    rdDepictor.Compute2DCoords(candidate)
                    smiles = Chem.MolToSmiles(candidate)
                    
                    if smiles and smiles not in successful_results:
                        successful_results.append(smiles)
                        print(f"✓✓✓ Attempt {trial} SUCCEEDED: {smiles[:80]}...")
                        
                except Exception as e:
                    print(f"✗✗✗ Attempt {trial} FAILED: {str(e)}")
                    import traceback
                    traceback.print_exc()
                    continue
            
            print(f"\n=== Generated {len(successful_results)} unique molecules ===\n")
            return successful_results if successful_results else [parent_smiles]
            
        except Exception as e:
            print(f"Error in replacement: {str(e)}")
            import traceback
            traceback.print_exc()
            return [parent_smiles]
    
    def _is_connected(self, mol: Chem.Mol) -> bool:
        """Check if molecule is fully connected (single fragment)."""
        return len(Chem.GetMolFrags(mol)) == 1
    
    def _score_alignment(self, mol_a: Chem.Mol, mol_b: Chem.Mol, atom_map: List[Tuple[int, int]]) -> float:
        """
        Helper to evaluate 2D overlap between molecules.
        
        Args:
            mol_a: First molecule
            mol_b: Second molecule
            atom_map: List of (atom_idx_a, atom_idx_b) tuples
            
        Returns:
            RMSD score (lower is better)
        """
        try:
            return rdMolAlign.AlignMol(mol_a, mol_b, atomMap=atom_map)
        except Exception:
            return float("inf")
    
    def generate_molecule_svg(self, smiles: str, width: int = 200, height: int = 200) -> str:
        """
        Generate an SVG visualization of a molecule from SMILES.
        
        Args:
            smiles: SMILES string representation of the molecule
            width: Width of the SVG in pixels
            height: Height of the SVG in pixels
            
        Returns:
            SVG string
            
        Raises:
            ValueError: If the SMILES string is invalid
        """
        mol = Chem.MolFromSmiles(smiles)
        
        if mol is None:
            raise ValueError(f"Invalid SMILES string: {smiles}")
        
        # Generate 2D coordinates
        AllChem.Compute2DCoords(mol)
        
        # Create the drawer
        drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
        
        # Set transparent background
        drawer.drawOptions().clearBackground = False
        
        # Draw molecule
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        
        svg = drawer.GetDrawingText()
        
        return svg
    
    def replace_substructure_with_llm(
        self,
        original_smiles: str,
        source_fragment_smiles: str,
        replacement_fragment_smiles: str
    ) -> Dict[str, Any]:
        """
        Use GPT to intelligently replace a substructure in a molecule.
        
        This function sends the original molecule, the fragment to replace,
        and the replacement fragment to GPT, which returns the resulting
        molecule SMILES with the replacement made.
        
        Args:
            original_smiles: SMILES of the original molecule
            source_fragment_smiles: SMILES of the fragment to be replaced
            replacement_fragment_smiles: SMILES of the replacement fragment
            
        Returns:
            Dictionary containing:
            - result_smiles: The resulting molecule SMILES
            - explanation: LLM explanation of what was done
            - success: Whether the replacement was successful
        """
        print(f"\n=== LLM Substructure Replacement ===")
        print(f"Original: {original_smiles}")
        print(f"Source fragment: {source_fragment_smiles}")
        print(f"Replacement fragment: {replacement_fragment_smiles}")
        
        # Validate inputs
        original_mol = Chem.MolFromSmiles(original_smiles)
        source_mol = Chem.MolFromSmiles(source_fragment_smiles)
        replacement_mol = Chem.MolFromSmiles(replacement_fragment_smiles)
        
        if not original_mol:
            raise ValueError(f"Invalid original SMILES: {original_smiles}")
        if not source_mol:
            raise ValueError(f"Invalid source fragment SMILES: {source_fragment_smiles}")
        if not replacement_mol:
            raise ValueError(f"Invalid replacement fragment SMILES: {replacement_fragment_smiles}")
        
        # Canonicalize inputs to avoid ambiguity for the LLM
        canonical_original = Chem.MolToSmiles(original_mol)
        canonical_source = Chem.MolToSmiles(source_mol)
        canonical_replacement = Chem.MolToSmiles(replacement_mol)
        print(f"Canonical Original: {canonical_original}")
        print(f"Canonical Source: {canonical_source}")
        print(f"Canonical Replacement: {canonical_replacement}")
        
        # Initialize LM Studio client (OpenAI-compatible API)
        # LM Studio runs locally on port 1234 by default
        lm_studio_url = os.getenv("LM_STUDIO_BASE_URL", "http://localhost:1234/v1")
        client = OpenAI(
            base_url=lm_studio_url,
            api_key="lm-studio"  # LM Studio doesn't require a real API key
        )
        
        # Construct the prompt for JSON structured output
        system_prompt = """You are an expert computational chemist specializing in molecular structure manipulation and SMILES notation.

Your task: Replace a specific substructure (source fragment) in a molecule with a different fragment (replacement fragment) while maintaining proper chemical connectivity and valence.

CRITICAL RULES:
1. The source fragment must be found and replaced in the original molecule
2. Maintain all chemical bonds and connectivity outside the replaced region
3. Ensure the replacement fragment connects at analogous positions to the source fragment
4. The resulting molecule must be chemically valid with correct valence for all atoms

You must respond with valid JSON in this exact format:
{
  "result_smiles": "the resulting molecule SMILES string after replacement"
}

Only output valid JSON. Do not include any explanation or additional text outside the JSON."""

        user_prompt = f"""Replace the source fragment with the replacement fragment in the original molecule.

Original Molecule SMILES: {canonical_original}
Source Fragment to Replace: {canonical_source}  
Replacement Fragment: {canonical_replacement}

Output your answer as JSON with a single field "result_smiles" containing the resulting molecule SMILES."""

        try:
            # Call local LM Studio
            # Note: LM Studio doesn't support response_format like OpenAI
            # We rely on prompt engineering to get JSON output
            print(f"Calling LM Studio at {lm_studio_url}...")
            response = client.chat.completions.create(
                model="gpt-oss-20b",  # Or whatever model name is loaded in LM Studio
                messages=[
                    {"role": "system", "content": system_prompt},
                    {"role": "user", "content": user_prompt}
                ],
                temperature=0.1,
                max_tokens=120000
            )
            
            response_text = response.choices[0].message.content.strip()
            print(f"LM Studio Response: {response_text}")
            
            # Parse JSON response
            result_smiles = None
            try:
                import json
                response_data = json.loads(response_text)
                result_smiles = response_data.get("result_smiles", "").strip()
                print(f"Parsed result_smiles from JSON: {result_smiles}")
            except json.JSONDecodeError as json_err:
                print(f"JSON parsing error: {json_err}")
                print(f"Raw response: {response_text}")
                
                # Fallback: try to extract SMILES from malformed response
                # Try "SMILES:" format
                if "SMILES:" in response_text:
                    smiles_match = re.search(r'SMILES:\s*([^\n,}]+)', response_text)
                    if smiles_match:
                        result_smiles = smiles_match.group(1).strip().strip('"')
                
                # Try JSON field extraction
                if not result_smiles:
                    smiles_match = re.search(r'"result_smiles"\s*:\s*"([^"]+)"', response_text)
                    if smiles_match:
                        result_smiles = smiles_match.group(1).strip()
                
                # Last resort: treat entire response as SMILES
                if not result_smiles:
                    candidate = response_text.strip().strip('`{}').strip('"').splitlines()[0].strip()
                    if candidate and len(candidate) > 10:
                        candidate_mol = Chem.MolFromSmiles(candidate, sanitize=False)
                        if candidate_mol:
                            print("Falling back to direct SMILES extraction.")
                            result_smiles = candidate
            
            print(f"Final extracted SMILES: {result_smiles}")
            
            if not result_smiles:
                print("SMILESSS:", result_smiles)
                raise ValueError(f"Could not extract SMILES from response: {response_text}")
            
            # Validate the result SMILES
            result_mol = Chem.MolFromSmiles(result_smiles, sanitize=False)
            if not result_mol:
                print(f"Warning: GPT returned invalid SMILES: {result_smiles}")
                raise ValueError(f"GPT returned invalid SMILES: {result_smiles}")
            
            # Canonicalize the result
            # canonical_smiles = Chem.MolToSmiles(result_mol)
            
            print(f"Result SMILES: {result_smiles}")
            
            return {
                "result_smiles": result_smiles,
                "explanation": None,
                "success": True
            }
            
        except Exception as e:
            print(f"Error in LLM replacement: {str(e)}")
            return {
                "result_smiles": original_smiles,
                "explanation": f"Error: {str(e)}",
                "success": False
            }
    
    def __del__(self):
        """Cleanup when controller is destroyed."""
        if hasattr(self, 'smart_chemist'):
            del self.smart_chemist

