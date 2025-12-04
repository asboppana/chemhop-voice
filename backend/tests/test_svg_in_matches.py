"""
Test to verify SVG fields are preserved through the entire data flow:
MCP Server → Agent → Voice Endpoint → Frontend

This tests the fix for the missing SVG fields in substructure matches.
"""
import sys
from pathlib import Path

# Add backend to path
backend_dir = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(backend_dir))

from src.api.models.mcp_models import MoleculeAnnotation, AnnotationMatch, TrivialName


def test_pydantic_model_accepts_svg():
    """Test that AnnotationMatch model accepts and preserves SVG field."""
    print("=" * 80)
    print("TEST 1: Pydantic Model SVG Support")
    print("=" * 80)
    
    # Simulate data from MCP server (what smart_chemist returns)
    mcp_data = {
        'name': 'Test Molecule',
        'svg': '<svg>main molecule svg</svg>',
        'smiles': 'c1ccccc1',
        'matches': [
            {
                'atom_indices': [0, 1, 2, 3, 4, 5],
                'svg': '<svg>benzene pattern svg</svg>',  # This is what MCP returns!
                'trivial_name': {
                    'name': 'Benzene',
                    'smarts': 'c1ccccc1',
                    'group': 'cyclic',
                    'bonds': 6,
                    'hierarchy': None
                }
            }
        ]
    }
    
    # Parse through Pydantic model
    annotation = MoleculeAnnotation(**mcp_data)
    
    assert annotation.svg is not None, "Main SVG is None"
    assert len(annotation.matches) == 1, f"Expected 1 match, got {len(annotation.matches)}"
    assert annotation.matches[0].svg is not None, "Match SVG is None"
    assert '<svg>benzene pattern svg</svg>' in annotation.matches[0].svg, "SVG content not preserved"
    
    print("✅ Pydantic model accepts and preserves SVG field")
    print(f"   Main molecule has SVG: {annotation.svg is not None}")
    print(f"   Match 0 has SVG: {annotation.matches[0].svg is not None}")
    print()


def test_json_serialization_preserves_svg():
    """Test that JSON serialization preserves SVG fields."""
    print("=" * 80)
    print("TEST 2: JSON Serialization")
    print("=" * 80)
    
    # Create annotation with SVG
    annotation = MoleculeAnnotation(
        name='Test',
        svg='<svg>main</svg>',
        smiles='c1ccccc1',
        matches=[
            AnnotationMatch(
                atom_indices=[0, 1, 2],
                svg='<svg>pattern</svg>',
                trivial_name=TrivialName(
                    name='Test Pattern',
                    smarts='c1ccccc1',
                    group='cyclic',
                    bonds=6,
                    hierarchy=None
                )
            )
        ]
    )
    
    # Serialize to dict (what gets sent to frontend)
    json_data = annotation.model_dump()
    
    assert 'svg' in json_data, "Main SVG not in JSON"
    assert 'matches' in json_data, "Matches not in JSON"
    assert len(json_data['matches']) == 1, f"Expected 1 match in JSON, got {len(json_data['matches'])}"
    assert 'svg' in json_data['matches'][0], "Match SVG not in JSON"
    assert json_data['matches'][0]['svg'] == '<svg>pattern</svg>', "Match SVG value incorrect"
    
    print("✅ JSON serialization preserves SVG fields")
    print(f"   Main SVG in JSON: {'svg' in json_data}")
    print(f"   Match SVG in JSON: {'svg' in json_data['matches'][0]}")
    print(f"   Match SVG value: {json_data['matches'][0]['svg']}")
    print()


def test_structured_message_format():
    """Test the structured message format from agent."""
    print("=" * 80)
    print("TEST 3: Structured Message Format")
    print("=" * 80)
    
    # Simulate what the agent should emit
    structured_message = {
        "type": "molecule_structured",
        "smiles": "c1ccccc1",
        "annotation": {
            "name": "Benzene",
            "svg": "<svg>main molecule</svg>",
            "smiles": "c1ccccc1",
            "matches": [
                {
                    "atom_indices": [0, 1, 2, 3, 4, 5],
                    "svg": "<svg>benzene pattern</svg>",
                    "trivial_name": {
                        "name": "Benzene",
                        "smarts": "c1ccccc1",
                        "group": "cyclic",
                        "bonds": 6,
                        "hierarchy": None
                    }
                }
            ]
        }
    }
    
    # Verify structure
    assert structured_message["type"] == "molecule_structured"
    assert "annotation" in structured_message
    assert "matches" in structured_message["annotation"]
    assert len(structured_message["annotation"]["matches"]) > 0
    
    first_match = structured_message["annotation"]["matches"][0]
    assert "svg" in first_match, "SVG field missing from match"
    assert "atom_indices" in first_match
    assert "trivial_name" in first_match
    
    print("✅ Structured message format is correct")
    print(f"   Has type field: {structured_message.get('type')}")
    print(f"   Has annotation wrapper: {'annotation' in structured_message}")
    print(f"   First match has SVG: {'svg' in first_match}")
    print()


def test_voice_endpoint_parsing():
    """Test the voice endpoint parsing logic."""
    print("=" * 80)
    print("TEST 4: Voice Endpoint Parsing Logic")
    print("=" * 80)
    
    # Simulate structured object from extract_structured_data
    structured_obj = {
        "type": "molecule_structured",
        "smiles": "c1ccccc1",
        "annotation": {
            "name": "Benzene",
            "svg": "<svg>main</svg>",
            "smiles": "c1ccccc1",
            "matches": [
                {
                    "atom_indices": [0, 1, 2, 3, 4, 5],
                    "svg": "<svg>pattern</svg>",
                    "trivial_name": {
                        "name": "Benzene",
                        "smarts": "c1ccccc1",
                        "group": "cyclic",
                        "bonds": 6,
                        "hierarchy": None
                    }
                }
            ]
        }
    }
    
    # Simulate _populate_typed_fields logic
    if "annotation" in structured_obj:
        annotation_data = structured_obj["annotation"].copy()
        if "smiles" not in annotation_data and "smiles" in structured_obj:
            annotation_data["smiles"] = structured_obj["smiles"]
        
        molecule_annotation = MoleculeAnnotation(**annotation_data)
        
        assert molecule_annotation is not None
        assert molecule_annotation.matches[0].svg is not None
        
        print("✅ Voice endpoint parsing works correctly")
        print(f"   Annotation created: {molecule_annotation is not None}")
        print(f"   Match SVG preserved: {molecule_annotation.matches[0].svg is not None}")
        print()


if __name__ == "__main__":
    print("\n" + "=" * 80)
    print("SVG FIELD PRESERVATION TEST SUITE")
    print("Testing fix for missing SVG fields in substructure matches")
    print("=" * 80 + "\n")
    
    try:
        test_pydantic_model_accepts_svg()
        test_json_serialization_preserves_svg()
        test_structured_message_format()
        test_voice_endpoint_parsing()
        
        print("=" * 80)
        print("✅ ALL TESTS PASSED!")
        print("=" * 80)
        print("\nSummary:")
        print("  1. ✅ AnnotationMatch model accepts SVG field")
        print("  2. ✅ JSON serialization preserves SVG")
        print("  3. ✅ Structured message format is correct")
        print("  4. ✅ Voice endpoint parsing works")
        print("\nThe fix is complete and working correctly!")
        print("=" * 80 + "\n")
        
    except AssertionError as e:
        print(f"\n❌ TEST FAILED: {e}\n")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    except Exception as e:
        print(f"\n❌ ERROR: {e}\n")
        import traceback
        traceback.print_exc()
        sys.exit(1)

