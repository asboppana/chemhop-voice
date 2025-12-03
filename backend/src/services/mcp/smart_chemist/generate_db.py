import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem
from tqdm import tqdm

# Add the smart_chemist directory to sys.path to allow imports
smart_chemist_dir = Path(__file__).resolve().parent
if str(smart_chemist_dir) not in sys.path:
    sys.path.insert(0, str(smart_chemist_dir))

from models import AnnotatedPattern, Base, engine, get_session  # noqa: E402

# Get the directory where this script is located
SCRIPT_DIR = Path(__file__).parent
DB_DIR = SCRIPT_DIR / "db"
CSV_PATH = DB_DIR / "smarts_with_hierarchy.csv"

# Drop all existing tables before recreating them
Base.metadata.drop_all(engine)
Base.metadata.create_all(engine)

# Get a database session
session = get_session()


# Initial population function
def populate_initial_patterns():
    patterns = []
    # read the smarts_with_hierarchy.csv file and parse it into a list of AnnotatedPattern objects
    smarts_with_hierarchy = pd.read_csv(CSV_PATH, skiprows=1)

    for index, row in tqdm(
        smarts_with_hierarchy.iterrows(),
        total=len(smarts_with_hierarchy),
        desc="Processing patterns",
    ):
        if row["group"] == "cyclic":
            mol = Chem.MolFromSmiles(row["SMARTS"])
        else:
            mol = Chem.MolFromSmarts(row["SMARTS"])
        if mol is None:
            print([x for x in row])
            continue
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
        n_heavy_atoms = mol.GetNumHeavyAtoms()
        n_other = (
            n_heavy_atoms
            - n_nitrogen
            - n_sulfur
            - n_oxygen
            - n_carbon
            - n_halogens
            - n_phospor
        )
        number_rings = len(Chem.GetSymmSSSR(mol))

        pattern = AnnotatedPattern(
            smarts=row["SMARTS"],
            trivial_name=row["trivialname"],
            group=row["group"],
            hierarchy=row["Hierarchy"],
            index_file=index + 1,
            heavy_atoms=n_heavy_atoms,
            num_rings=number_rings,
            n_nitrogens=n_nitrogen,
            n_sulfur=n_sulfur,
            n_oxygen=n_oxygen,
            n_carbon=n_carbon,
            n_phosphor=n_phospor,
            n_halogens=n_halogens,
            n_other_atom=n_other,
        )
        
        patterns.append(pattern)

    session.add_all(patterns)
    session.commit()
    print("Database populated with initial AnnotatedPattern entries.")


if __name__ == "__main__":
    populate_initial_patterns()
