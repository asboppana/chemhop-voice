import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
import requests
import xml.etree.ElementTree as ET
import pandas as pd
from pathlib import Path

# Import the database models using relative imports
from ..models import AnnotatedPattern, get_session

def convert_string_input_to_smiles(input_string):
    """Parse an input request string."""
    if Chem.MolFromSmiles(input_string) is not None:
        return [input_string]
    elif input_string.startswith("chembl:"):
        # example: chembl:CHEMBL50894
        chembl_id = input_string.removeprefix("chembl:")
        url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/{chembl_id}.json"
        res = requests.get(url).json()
        return [res["molecule_structures"]["canonical_smiles"]]
    elif input_string.startswith("chebi"):
        # example: chebi:138488
        chebi_id = input_string.removeprefix("chebi:")
        chebi_url = f"http://www.ebi.ac.uk/webservices/chebi/2.0/test/getCompleteEntity?chebiId={chebi_id}"
        response = requests.get(chebi_url)
        root = ET.fromstring(response.text)
        return [root.find(".//{*}smiles").text]
    elif input_string.startswith("pubchem"):
        # example: pubchem:5005498
        cid = input_string.removeprefix("pubchem:")
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/JSON"
        data = requests.get(url).json()
        for property in data["PC_Compounds"][0]["props"]:
            if property["urn"]["label"] == "SMILES" and property["urn"]["name"] == "Absolute":
                return [property["value"]["sval"]]
        return [""]
    else:
        # Query the database for matching patterns
        session = get_session()
        try:
            patterns = session.query(AnnotatedPattern).filter(
                AnnotatedPattern.trivial_name == input_string
            ).all()
            
            for pattern in patterns:
                if pattern.group == "cyclic":
                    return [f"{pattern.smarts} {pattern.trivial_name}", "CCCCCCC"]
                else:
                    # Try to find test cases if available
                    testcases_path = Path(__file__).parent.parent / "db" / "smarts_testcases.csv"
                    if testcases_path.exists():
                        testdata = pd.read_csv(testcases_path)
                        matching_data = testdata.loc[testdata["pattern"] == pattern.trivial_name]
                        if matching_data.shape[0] == 0:
                            return [""]
                        else:
                            molecule_list = []
                            for index, row in matching_data.iterrows():
                                molecule_list.append(f"{row['smiles']} {row['comment']}")
                            return molecule_list
                    return [""]
            return [""]
        finally:
            session.close()


def check_subset_of_lists(list1, list2):
    """Check if list1 is a proper subset of list2."""
    if len(list1) >= len(list2):
        return False
    set1 = set(list1)
    set2 = set(list2)
    return set1.issubset(set2)

def check_lists_equal(list1, list2):
    """Check if list1 and list2 contain the same elements."""
    if len(list1) != len(list2):
        return False
    return set(list1) == set(list2)

def check_list_in_list(list1, list2):
    """Check if all elements in list1 are present in list2."""
    set2 = set(list2)
    return all(element in set2 for element in list1)


class SmartChemist:
    pattern_dictionary = {}  # Dictionary with index as key and built mol object as value
    # This dictionary should drastically reduce the time needed for multi-mol queries
    
    def __init__(self):
        """Initialize SmartChemist with a database session."""
        self.session = get_session()
    
    def __del__(self):
        """Close the database session when the object is destroyed."""
        if hasattr(self, 'session'):
            self.session.close()
    
    @staticmethod
    def _check_hierarchy(matches: list):
        for match in matches:
            hierarchy = match["trivial_name"]["hierarchy"]
            atoms_match = match["atom_indices"]
            if hierarchy and hierarchy != "[]":
                try:
                    hierarchy_pattern_indexes = [int(x) for x in hierarchy.replace("[", "").replace("]", "").split(",") if x.strip()]
                    for submatch in matches:
                        atoms_submatch = submatch["atom_indices"]
                        submatch_id = submatch["trivial_name"].get("index")
                        if submatch_id is None or submatch_id == "":
                            continue
                        if submatch_id in hierarchy_pattern_indexes and atoms_submatch[0] in atoms_match:
                            submatch["trivial_name"]["group"] = "overshadowed"
                except (ValueError, AttributeError):
                    # Skip if hierarchy parsing fails
                    pass

    @staticmethod
    def _check_overshadowed_patterns(matches: list):
        SmartChemist._check_hierarchy(matches)
        for match in matches:
            if match["trivial_name"]["group"] == "overshadowed":
                continue
            for submatch in matches:
                atoms_match = match["atom_indices"]
                atoms_submatch = submatch["atom_indices"]
                if check_subset_of_lists(atoms_match, atoms_submatch):
                    match["trivial_name"]["group"] = "overshadowed"
                    break
                elif check_lists_equal(atoms_match, atoms_submatch):
                    bonds1 = match["trivial_name"]["bonds"]
                    bonds2 = submatch["trivial_name"]["bonds"]
                    if bonds1 > bonds2:
                        submatch["trivial_name"]["group"] = "overshadowed"
                    elif bonds2 < bonds1:
                        match["trivial_name"]["group"] = "overshadowed"

    def _match_smarts_patterns(self, mol: rdkit.Chem.rdchem.Mol, remove_overshadowed_patterns: bool = False):
        matches = []
        # iterate all annotated SMARTS patterns from our database
        heavy_atoms = mol.GetNumHeavyAtoms()
        number_rings = len(Chem.GetSymmSSSR(mol))
        n_nitrogen = n_sulfur = n_oxygen = n_carbon = n_halogens = n_phospor = 0
        for atom in mol.GetAtoms():
            number = atom.GetAtomicNum()
            if number == 6:
                n_carbon += 1
            elif number == 7:
                n_nitrogen += 1
            elif number == 8:
                n_oxygen += 1
            elif number == 15:
                n_phospor += 1
            elif number == 16:
                n_sulfur += 1
            elif number == 9 or number == 17 or number == 35 or number == 53:
                n_halogens += 1
        n_other = heavy_atoms - n_nitrogen - n_sulfur - n_oxygen - n_carbon - n_halogens - n_phospor
        
        # Query database with filters to improve performance
        db_patterns = self.session.query(AnnotatedPattern).filter(
            AnnotatedPattern.heavy_atoms <= heavy_atoms,
            AnnotatedPattern.num_rings <= number_rings,
            AnnotatedPattern.n_nitrogens <= n_nitrogen,
            AnnotatedPattern.n_oxygen <= n_oxygen,
            AnnotatedPattern.n_sulfur <= n_sulfur,
            AnnotatedPattern.n_carbon <= n_carbon,
            AnnotatedPattern.n_halogens <= n_halogens,
            AnnotatedPattern.n_phosphor <= n_phospor,
            AnnotatedPattern.n_other_atom <= n_other
        ).all()
        
        for db_row in db_patterns:
            if db_row.id in self.pattern_dictionary:
                pattern = self.pattern_dictionary[db_row.id]
            else:
                if db_row.group == "cyclic":
                    pattern = Chem.MolFromSmiles(db_row.smarts)
                else:
                    pattern = Chem.MolFromSmarts(db_row.smarts)
                self.pattern_dictionary[db_row.id] = pattern
            
            if pattern is None:
                continue
                
            if mol.HasSubstructMatch(pattern, useChirality=True):
                hit_atom_indices_list = mol.GetSubstructMatches(pattern, useChirality=True)
                for hit_atom_indices in hit_atom_indices_list:
                    matches.append(
                        {
                            "atom_indices": hit_atom_indices,
                            "trivial_name": {
                                "name": db_row.trivial_name,
                                "smarts": db_row.smarts,
                                "group": db_row.group,
                                "bonds": pattern.GetNumBonds(),
                                "hierarchy": db_row.hierarchy,
                                "index": db_row.index_file,
                            },
                        }
                    )
        self._check_overshadowed_patterns(matches)
        if remove_overshadowed_patterns:
            matches = [x for x in matches if x["trivial_name"]["group"] != "overshadowed"]
        return matches

    @staticmethod
    def _mol_to_image_str(mol: rdkit.Chem.rdchem.Mol, width: int, height: int) -> str:
        AllChem.Compute2DCoords(mol)
        if not mol.GetConformer(0):
            raise ValueError("Failed to generate 2D conformation for mol")
        drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        result = drawer.GetDrawingText()
        return result

    def mol_to_annotation_json(self, mol: rdkit.Chem.rdchem.Mol) -> dict:
        # search the database for annotations of substructures
        db_matches = self._match_smarts_patterns(mol)
        if mol.HasProp("_Name"):
            name = mol.GetProp("_Name")
        else:
            name = "No Name"
        return {
            "name": name,
            "svg": self._mol_to_image_str(mol, 400, 400),
            "matches": db_matches,
            "smiles": Chem.MolToSmiles(mol),
        }
