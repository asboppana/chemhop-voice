"""
Controller for molecule analysis operations.

Handles the business logic for analyzing molecules using SMILES strings
and interfacing with the SmartChemist tool.
"""
from typing import List
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import rdMolDraw2D

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
    
    def __del__(self):
        """Cleanup when controller is destroyed."""
        if hasattr(self, 'smart_chemist'):
            del self.smart_chemist

