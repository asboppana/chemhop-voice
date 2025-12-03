# Smart Chemist - SMARTS Pattern Matching

A chemistry tool for annotating molecules using SMARTS patterns and functional group recognition.

## Overview

Smart Chemist analyzes molecular structures (provided as SMILES strings) and identifies functional groups and structural patterns using a database of annotated SMARTS patterns.

## Project Structure

```
smart_chemist/
├── README.md                  # This file
├── __init__.py               # Package initialization
├── db/                       # Database files
│   ├── smart_chemist.db      # SQLite database with patterns
│   └── smarts_with_hierarchy.csv  # Source data for patterns
├── models/                   # SQLAlchemy models
│   └── __init__.py          # AnnotatedPattern model and database setup
├── tools/                    # Main tools
│   ├── __init__.py
│   └── smart_chemist.py     # SmartChemist class for pattern matching
├── generate_db.py           # Script to generate/regenerate the database
└── test_query.py            # Test script to verify functionality
```

## Database Setup

### Initial Setup

To create the database from the CSV file:

```bash
cd /path/to/backend/src/services/mcp/smart_chemist
python -m smart_chemist.generate_db
```

This will:
1. Drop any existing tables (clean slate)
2. Create the `annotated_pattern` table
3. Parse the CSV file and populate the database with ~40,000+ patterns
4. Calculate molecular descriptors for each pattern

### Database Schema

The `AnnotatedPattern` table contains:
- `id`: Primary key
- `trivial_name`: Name of the pattern (e.g., "Benzene", "Carboxylic acid")
- `smarts`: SMARTS pattern string
- `group`: Pattern group (e.g., "cyclic", "functional_group")
- `index_file`: Index from source file
- `heavy_atoms`: Number of heavy atoms
- `num_rings`: Number of rings
- `n_nitrogens`, `n_oxygen`, `n_sulfur`, `n_carbon`, `n_halogens`, `n_phosphor`, `n_other_atom`: Element counts
- `hierarchy`: Hierarchical relationship to other patterns

## Usage

### Basic Usage

```python
from src.services.mcp import SmartChemist
from rdkit import Chem

# Initialize the chemist
chemist = SmartChemist()

# Analyze a molecule (aspirin example)
mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
result = chemist.mol_to_annotation_json(mol)

print(f"Found {len(result['matches'])} matches")
for match in result['matches']:
    print(f"- {match['trivial_name']['name']}")
```

### Querying the Database Directly

```python
from src.services.mcp import AnnotatedPattern, get_session

session = get_session()

# Query patterns with rings
patterns = session.query(AnnotatedPattern).filter(
    AnnotatedPattern.num_rings >= 1
).limit(10).all()

for pattern in patterns:
    print(f"{pattern.trivial_name}: {pattern.smarts}")

session.close()
```

### Converting Input Strings

The `convert_string_input_to_smiles` function supports multiple input formats:

```python
from src.services.mcp.smart_chemist.tools.smart_chemist import convert_string_input_to_smiles

# SMILES string
convert_string_input_to_smiles("CC(=O)O")  # Returns: ["CC(=O)O"]

# ChEMBL ID
convert_string_input_to_smiles("chembl:CHEMBL50894")

# ChEBI ID
convert_string_input_to_smiles("chebi:138488")

# PubChem CID
convert_string_input_to_smiles("pubchem:5005498")

# Pattern name from database
convert_string_input_to_smiles("Benzene")
```

## Testing

Run the test suite:

```bash
cd /path/to/backend/src/services/mcp/smart_chemist
python -m smart_chemist.test_query
```

This will:
1. Count total patterns in database
2. Query sample patterns
3. Test the SmartChemist class with aspirin
4. Display matches found

## Features

### Pattern Matching
- Matches molecular substructures using SMARTS patterns
- Supports chirality in matching
- Identifies all instances of patterns in a molecule

### Hierarchy & Overshadowing
- Patterns can have hierarchical relationships
- More specific patterns can overshadow general ones
- Equal patterns are compared by number of bonds

### Performance Optimization
- Pattern dictionary caches compiled molecular patterns
- Database queries filter by molecular descriptors before matching
- Only relevant patterns are tested against each molecule

## Dependencies

- `sqlalchemy>=2.0.44`: Database ORM
- `rdkit>=2025.9.3`: Chemistry toolkit
- `pandas>=2.3.3`: CSV processing
- `tqdm>=4.67.1`: Progress bars
- `requests>=2.32.5`: External API calls

## Notes

- The database file (`smart_chemist.db`) is regenerated from scratch each time `generate_db.py` is run
- Duplicate pattern names in the source CSV are allowed (unique constraint removed)
- The database uses SQLite for simplicity and portability
- All paths are resolved relative to the module location for flexibility

