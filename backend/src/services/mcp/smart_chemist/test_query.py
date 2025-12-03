"""Test script to verify SQLAlchemy queries work correctly."""

import sys
from pathlib import Path

# Add the smart_chemist directory to sys.path to allow imports
smart_chemist_dir = Path(__file__).resolve().parent
if str(smart_chemist_dir) not in sys.path:
    sys.path.insert(0, str(smart_chemist_dir))

from models import AnnotatedPattern, get_session
from rdkit import Chem

def test_database_query():
    """Test basic database queries."""
    session = get_session()
    
    try:
        # Test 1: Count all patterns
        total_patterns = session.query(AnnotatedPattern).count()
        print(f"✓ Total patterns in database: {total_patterns}")
        
        # Test 2: Get first 5 patterns
        first_patterns = session.query(AnnotatedPattern).limit(5).all()
        print("\n✓ First 5 patterns:")
        for pattern in first_patterns:
            print(f"  - ID: {pattern.id}, Name: {pattern.trivial_name}, Group: {pattern.group}")
        
        # Test 3: Query by specific criteria
        aromatic_patterns = session.query(AnnotatedPattern).filter(
            AnnotatedPattern.num_rings >= 1
        ).limit(5).all()
        print("\n✓ First 5 patterns with rings:")
        for pattern in aromatic_patterns:
            print(f"  - Name: {pattern.trivial_name}, Rings: {pattern.num_rings}, Heavy Atoms: {pattern.heavy_atoms}")
        
        # Test 4: Test SmartChemist class
        print("\n✓ Testing SmartChemist class...")
        
        from tools.smart_chemist import SmartChemist
        import json
        
        chemist = SmartChemist()
        
        # Test with a more complex molecule containing multiple functional groups
        # Aspirin: contains aromatic ring, carboxylic acid, and ester groups
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        if mol:
            result = chemist.mol_to_annotation_json(mol)
            print( " - Analyzed Aspirin (CC(=O)Oc1ccccc1C(=O)O)")
            print(f" - Found {len(result['matches'])} matches")
            if result['matches']:
                print(f"  - First match: {result['matches'][0]['trivial_name']['name']}")
            
            print(json.dumps(result, indent=4))    
            
        print("\n✅ All tests passed!")
        
    except Exception as e:
        print(f"\n❌ Error: {e}")
        import traceback
        traceback.print_exc()
    finally:
        session.close()


if __name__ == "__main__":
    test_database_query()

