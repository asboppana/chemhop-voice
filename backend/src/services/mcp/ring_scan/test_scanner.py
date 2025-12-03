#!/usr/bin/env python3
"""
Test script for bio-isostere scanner.

Tests the bio-isostere scanning functionality with various ring structures.
"""
import sys
from pathlib import Path

# Add backend to path for imports
backend_path = Path(__file__).parent.parent.parent.parent
sys.path.insert(0, str(backend_path))

from src.services.mcp.ring_scan import get_scanner


def test_scanner():
    """Test the bio-isostere scanner with example rings."""
    
    print("=" * 80)
    print("Bio-isostere Scanner Test Suite")
    print("=" * 80)
    print()
    
    # Initialize scanner
    print("Initializing scanner...")
    try:
        scanner = get_scanner()
        print(f"✓ Scanner initialized successfully")
        print(f"  - ertl_data: {len(scanner.ertl_data)} molecules")
        print(f"  - chemspace: {len(scanner.chemspace_data)} molecules")
        print(f"  - clusters: {len(scanner.centroids)} centroids")
        print(f"  - ring_smiles: {len(scanner.ring_smiles)} rings")
        print()
    except Exception as e:
        print(f"✗ Failed to initialize scanner: {e}")
        return False
    
    # Test cases: various ring structures
    test_cases = [
        ("Benzene", "c1ccccc1", "6-membered aromatic"),
        ("Naphthalene", "c1ccc2ccccc2c1", "Fused bicyclic aromatic"),
        ("Pyridine", "c1ccncc1", "6-membered N-heterocycle"),
        ("Thiazole", "c1scnc1", "5-membered S,N-heterocycle"),
        ("Furan", "c1ccoc1", "5-membered O-heterocycle"),
        ("Imidazole", "c1[nH]cnc1", "5-membered N,N-heterocycle"),
        ("Cyclohexane", "C1CCCCC1", "6-membered saturated"),
        ("Piperidine", "C1CCNCC1", "6-membered saturated N-heterocycle"),
        ("Quinoline", "c1ccc2ncccc2c1", "Fused bicyclic N-heterocycle"),
        ("Pyrimidine", "c1cncnc1", "6-membered N,N-heterocycle"),
    ]
    
    all_passed = True
    
    for name, smiles, description in test_cases:
        print("-" * 80)
        print(f"Test: {name} ({smiles})")
        print(f"Description: {description}")
        print("-" * 80)
        
        try:
            # Scan for bio-isosteres with adaptive scoring
            results = scanner.scan(
                query_smiles=smiles,
                top_k=10,
                min_similarity=0.2,
                use_bio_filters=True
            )
            
            if len(results) == 0:
                print(f"⚠ Warning: No matches found for {name}")
            else:
                print(f"✓ Found {len(results)} bio-isosteres")
                
                # Count by source
                sources = {}
                for r in results:
                    source = r.get('source', 'unknown')
                    sources[source] = sources.get(source, 0) + 1
                
                print(f"  Sources: {', '.join(f'{k}={v}' for k, v in sources.items())}")
            
            # Display top 5 matches
            print(f"\n  Top 5 Matches:")
            for i, match in enumerate(results[:5], 1):
                source = match.get('source', '?')
                bio_score = match.get('bio_isostere_score', 0)
                morgan_sim = match.get('similarity', 0)
                pharm_sim = match.get('pharmacophore_similarity', 0)
                topo_sim = match.get('topology_similarity', 0)
                centroid_smiles = match.get('centroid_smiles', 'N/A')
                
                print(f"\n  {i}. [{source.upper()}] {centroid_smiles}")
                print(f"     Bio-Score: {bio_score:.3f} (M={morgan_sim:.3f}, P={pharm_sim:.3f}, T={topo_sim:.3f})")
                
                # Show property differences if available
                delta_props = match.get('delta_properties', {})
                if delta_props:
                    delta_logp = delta_props.get('delta_logp', 0)
                    delta_tpsa = delta_props.get('delta_tpsa', 0)
                    print(f"     Δ logP: {delta_logp:.2f}, Δ TPSA: {delta_tpsa:.1f} Ų")
            
            print()
            
        except Exception as e:
            print(f"✗ Error scanning {name}: {e}")
            import traceback
            traceback.print_exc()
            all_passed = False
            print()
    
    # Test cluster info retrieval - collect unique cluster IDs from scan results
    print("=" * 80)
    print("Testing Cluster Info Retrieval")
    print("=" * 80)
    
    # Collect cluster IDs found during scans (from cluster source matches)
    cluster_ids_to_test = set()
    for name, smiles, _ in test_cases[:5]:  # Test first 5
        try:
            results = scanner.scan(query_smiles=smiles, top_k=10, min_similarity=0.2, use_bio_filters=True)
            for r in results:
                if r.get('source') == 'clusters' and r.get('cluster_id') is not None:
                    cluster_ids_to_test.add(r['cluster_id'])
                    if len(cluster_ids_to_test) >= 3:  # Test up to 3 different clusters
                        break
            if len(cluster_ids_to_test) >= 3:
                break
        except Exception:
            pass
    
    # Fallback to cluster 0 if none found
    if not cluster_ids_to_test:
        cluster_ids_to_test = {0}
    
    print(f"Testing {len(cluster_ids_to_test)} cluster(s): {sorted(cluster_ids_to_test)}")
    print()
    
    for cluster_id in sorted(cluster_ids_to_test):
        try:
            info = scanner.get_cluster_info(cluster_id)
            if info:
                print(f"  ✓ Cluster {cluster_id}:")
                print(f"    - Members: {info['num_members']}")
                if info['member_smiles']:
                    print(f"    - First 3: {', '.join(info['member_smiles'][:3])}")
            else:
                print(f"  ✗ Cluster {cluster_id} not found")
                all_passed = False
        except Exception as e:
            print(f"  ✗ Error retrieving cluster {cluster_id}: {e}")
            all_passed = False
    
    print()
    print("=" * 80)
    if all_passed:
        print("✅ All tests PASSED!")
    else:
        print("⚠ Some tests had warnings or errors")
    print("=" * 80)
    
    return all_passed


if __name__ == "__main__":
    try:
        success = test_scanner()
        sys.exit(0 if success else 1)
    except KeyboardInterrupt:
        print("\n\nTest interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\n\nUnexpected error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
