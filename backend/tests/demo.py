"""
Comprehensive demo for ChatAgent MCP integration.

Tests the full scaffold-hopping workflow:
1. Identifier -> SMILES conversion -> Annotation
2. Multi-turn conversations with group selection
3. Bioisostere scanning with AND/OR logic
4. Edge cases and error handling

Prerequisites:
1. Run MCP server: cd backend && python mcp_start.py
2. Set MCP_SERVER_URL if using ngrok
3. Run: OPENAI_API_KEY=your_key CHEM_AGENT_DEBUG=1 python tests/demo.py
"""

import sys
from pathlib import Path

# Add backend to path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.services.chatAgent.agent import ChatAgent


def print_section(title: str):
    """Print a formatted section header."""
    print("\n" + "=" * 80)
    print(f"📋 {title}")
    print("=" * 80)


def print_response(response: str, label: str = "Response"):
    """Print a formatted response."""
    print(f"\n✅ {label}:")
    print("-" * 80)
    print(response)
    print("-" * 80)


def scenario_1_simple_annotation():
    """Scenario 1: Simple molecule annotation (no follow-up)."""
    print_section("SCENARIO 1: Simple Molecule Annotation")
    
    agent = ChatAgent()
    
    print("\n👤 User: 'Analyze aspirin (CHEMBL25) and show me its functional groups'")
    response = agent.chat("Analyze aspirin (CHEMBL25) and show me its functional groups")
    print_response(response, "Agent")
    
    print("\n💡 Expected: Should convert CHEMBL25 -> SMILES -> annotate functional groups")


def scenario_2_full_workflow_with_single_group():
    """Scenario 2: Full workflow - annotation then modify a single group."""
    print_section("SCENARIO 2: Full Workflow - Modify Single Group")
    
    agent = ChatAgent()
    history = []
    
    # Step 1: Initial molecule
    print("\n👤 User: 'I have benzene c1ccccc1, can you annotate it?'")
    msg1 = "I have benzene c1ccccc1, can you annotate it?"
    resp1 = agent.chat(msg1, prior=history)
    history.extend([
        {"role": "user", "content": msg1},
        {"role": "assistant", "content": resp1}
    ])
    print_response(resp1, "Agent (Step 1)")
    
    # Step 2: Ask for bioisosteres for the benzene ring
    print("\n👤 User: 'Find bioisosteres for the benzene ring'")
    msg2 = "Find bioisosteres for the benzene ring"
    resp2 = agent.chat(msg2, prior=history)
    print_response(resp2, "Agent (Step 2)")
    
    print("\n💡 Expected: Should scan for benzene bioisosteres and return alternatives")


def scenario_3_multi_group_AND():
    """Scenario 3: Modify multiple groups with AND logic (Cartesian product)."""
    print_section("SCENARIO 3: Multi-Group Modification with AND Logic")
    
    agent = ChatAgent()
    history = []
    
    # Step 1: Analyze aspirin
    print("\n👤 User: 'Analyze aspirin CC(=O)Oc1ccccc1C(=O)O'")
    msg1 = "Analyze aspirin CC(=O)Oc1ccccc1C(=O)O"
    resp1 = agent.chat(msg1, prior=history)
    history.extend([
        {"role": "user", "content": msg1},
        {"role": "assistant", "content": resp1}
    ])
    print_response(resp1, "Agent (Step 1 - Annotation)")
    
    # Step 2: Request modifications with AND logic
    print("\n👤 User: 'Find bioisosteres for the acetyl group AND the benzene ring'")
    msg2 = "Find bioisosteres for the acetyl group AND the benzene ring"
    resp2 = agent.chat(msg2, prior=history)
    print_response(resp2, "Agent (Step 2 - AND Logic)")
    
    print("\n💡 Expected: Should find bioisosteres for BOTH groups")
    print("   - If acetyl has 3 alternatives and benzene has 5 alternatives")
    print("   - AND logic = 3 × 5 = 15 combinations (Cartesian product)")


def scenario_4_multi_group_OR():
    """Scenario 4: Modify multiple groups with OR logic (union)."""
    print_section("SCENARIO 4: Multi-Group Modification with OR Logic")
    
    agent = ChatAgent()
    history = []
    
    # Step 1: Analyze molecule
    print("\n👤 User: 'Annotate this molecule: CC(=O)Oc1ccccc1C(=O)O'")
    msg1 = "Annotate this molecule: CC(=O)Oc1ccccc1C(=O)O"
    resp1 = agent.chat(msg1, prior=history)
    history.extend([
        {"role": "user", "content": msg1},
        {"role": "assistant", "content": resp1}
    ])
    print_response(resp1, "Agent (Step 1)")
    
    # Step 2: Request modifications with OR logic
    print("\n👤 User: 'Find bioisosteres for the benzene ring OR the carboxylic acid group'")
    msg2 = "Find bioisosteres for the benzene ring OR the carboxylic acid group"
    resp2 = agent.chat(msg2, prior=history)
    print_response(resp2, "Agent (Step 2 - OR Logic)")
    
    print("\n💡 Expected: Should find bioisosteres for EITHER group")
    print("   - If benzene has 5 alternatives and carboxylic acid has 4 alternatives")
    print("   - OR logic = 5 + 4 = 9 total alternatives (union)")


def scenario_5_complex_boolean():
    """Scenario 5: Complex boolean logic - A OR (B AND C)."""
    print_section("SCENARIO 5: Complex Boolean Logic")
    
    agent = ChatAgent()
    history = []
    
    # Step 1: Analyze molecule
    print("\n👤 User: 'What functional groups are in aspirin?'")
    msg1 = "What functional groups are in aspirin?"
    resp1 = agent.chat(msg1, prior=history)
    history.extend([
        {"role": "user", "content": msg1},
        {"role": "assistant", "content": resp1}
    ])
    print_response(resp1, "Agent (Step 1)")
    
    # Step 2: Complex boolean
    print("\n👤 User: 'Find bioisosteres for the benzene ring OR (acetyl AND carboxylic acid)'")
    msg2 = "Find bioisosteres for the benzene ring OR (acetyl AND carboxylic acid)"
    resp2 = agent.chat(msg2, prior=history)
    print_response(resp2, "Agent (Step 2 - Complex Logic)")
    
    print("\n💡 Expected: benzene alternatives OR (acetyl × carboxylic acid combinations)")


def scenario_6_unclear_selection():
    """Scenario 6: User doesn't specify which group - agent should ask."""
    print_section("SCENARIO 6: Unclear Group Selection (Edge Case)")
    
    agent = ChatAgent()
    history = []
    
    # Step 1: Analyze molecule
    print("\n👤 User: 'Analyze CHEMBL25'")
    msg1 = "Analyze CHEMBL25"
    resp1 = agent.chat(msg1, prior=history)
    history.extend([
        {"role": "user", "content": msg1},
        {"role": "assistant", "content": resp1}
    ])
    print_response(resp1, "Agent (Step 1)")
    
    # Step 2: Vague request
    print("\n👤 User: 'Find some bioisosteres'")
    msg2 = "Find some bioisosteres"
    resp2 = agent.chat(msg2, prior=history)
    print_response(resp2, "Agent (Step 2)")
    
    print("\n💡 Expected: Agent should ask user to specify which group/section to modify")


def scenario_7_non_sequitur():
    """Scenario 7: User asks something unrelated after annotation."""
    print_section("SCENARIO 7: Non-Sequitur (Edge Case)")
    
    agent = ChatAgent()
    history = []
    
    # Step 1: Analyze molecule
    print("\n👤 User: 'Annotate benzene'")
    msg1 = "Annotate benzene"
    resp1 = agent.chat(msg1, prior=history)
    history.extend([
        {"role": "user", "content": msg1},
        {"role": "assistant", "content": resp1}
    ])
    print_response(resp1, "Agent (Step 1)")
    
    # Step 2: Unrelated question
    print("\n👤 User: 'What's the weather like?'")
    msg2 = "What's the weather like?"
    resp2 = agent.chat(msg2, prior=history)
    print_response(resp2, "Agent (Step 2)")
    
    print("\n💡 Expected: Agent should politely redirect to chemistry/drug discovery topics")


def scenario_8_name_to_smiles():
    """Scenario 8: Convert drug name to SMILES."""
    print_section("SCENARIO 8: Drug Name to SMILES Conversion")
    
    agent = ChatAgent()
    
    print("\n👤 User: 'What's the SMILES for ibuprofen?'")
    response = agent.chat("What's the SMILES for ibuprofen?")
    print_response(response, "Agent")
    
    print("\n💡 Expected: Should use get_smiles_from_name or convert_identifier_to_smiles")


def scenario_9_properties_request():
    """Scenario 9: Request ADMET properties for modified molecules."""
    print_section("SCENARIO 9: Request Properties for Candidates")
    
    agent = ChatAgent()
    history = []
    
    # Step 1: Get bioisosteres
    print("\n👤 User: 'Find bioisosteres for pyridine'")
    msg1 = "Find bioisosteres for pyridine"
    resp1 = agent.chat(msg1, prior=history)
    history.extend([
        {"role": "user", "content": msg1},
        {"role": "assistant", "content": resp1}
    ])
    print_response(resp1, "Agent (Step 1)")
    
    # Step 2: Request properties
    print("\n👤 User: 'What are the ADMET properties for the top 3 candidates?'")
    msg2 = "What are the ADMET properties for the top 3 candidates?"
    resp2 = agent.chat(msg2, prior=history)
    print_response(resp2, "Agent (Step 2)")
    
    print("\n💡 Expected: Should call predict_admet_properties for each candidate")


def main():
    """Run all demo scenarios."""
    print("\n" + "🧪 " * 20)
    print("ChatAgent MCP Comprehensive Demo")
    print("Testing Full Scaffold-Hopping Workflow")
    print("🧪 " * 20)
    
    scenarios = [
        ("Simple Cases", [
            scenario_1_simple_annotation,
            scenario_8_name_to_smiles,
        ]),
        ("Full Workflow", [
            scenario_2_full_workflow_with_single_group,
            scenario_3_multi_group_AND,
            scenario_4_multi_group_OR,
            scenario_5_complex_boolean,
        ]),
        ("Properties & ADMET", [
            scenario_9_properties_request,
        ]),
        ("Edge Cases", [
            scenario_6_unclear_selection,
            scenario_7_non_sequitur,
        ]),
    ]
    
    for category, funcs in scenarios:
        print(f"\n{'='*80}")
        print(f"📚 {category.upper()}")
        print(f"{'='*80}")
        
        for func in funcs:
            try:
                func()
            except Exception as e:
                print(f"\n❌ Error in {func.__name__}: {e}")
                import traceback
                traceback.print_exc()
            
            input("\n⏸️  Press Enter to continue to next scenario...")
    
    print("\n" + "=" * 80)
    print("✨ Demo Complete!")
    print("=" * 80)


if __name__ == "__main__":
    main()
