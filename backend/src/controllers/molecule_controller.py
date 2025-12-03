"""
Controller for molecule analysis operations.

Handles the business logic for analyzing molecules using SMILES strings
and interfacing with the SmartChemist tool.
"""
from rdkit import Chem

from src.services.mcp.smart_chemist.tools.smart_chemist import SmartChemist
from src.api.models.molecule import MoleculeResponse, Match, TrivialName


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
    
    def __del__(self):
        """Cleanup when controller is destroyed."""
        if hasattr(self, 'smart_chemist'):
            del self.smart_chemist

