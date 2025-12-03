# Bio-isostere Scanner

A powerful algorithm for finding bio-isosteric ring replacements by scanning input ring SMILES against ~34,000 cluster representatives from a database of ~4 million medicinal chemistry relevant ring systems.

## Overview

The bio-isostere scanner uses **Tanimoto similarity** with **Morgan fingerprints** to identify structurally similar ring systems that could serve as bio-isosteric replacements in drug design. The algorithm compares query rings against cluster centroids from a pre-clustered database, enabling fast similarity search across millions of ring structures.

## Database

- **Total rings**: 3,931,782 medicinal chemistry relevant ring systems
- **Cluster centroids**: 34,389 representative structures
- **Source**: Database of 4 million medicinal chemistry relevant ring systems by Peter Ertl (J.Chem.Inf.Model. 2024)
- **Clustering method**: BIRCH algorithm with Tanimoto similarity (threshold 0.3)

## API Endpoints

### 1. Scan for Bio-isosteres

**POST** `/molecule/bioisostere/scan`

Scan for bio-isosteric ring replacements.

#### Request Body

```json
{
  "ring_smiles": "c1ccccc1",
  "top_k": 20,
  "min_similarity": 0.3
}
```

**Parameters:**
- `ring_smiles` (string, required): SMILES string of the query ring structure
- `top_k` (integer, optional): Number of top matches to return (1-100, default: 20)
- `min_similarity` (float, optional): Minimum Tanimoto similarity threshold (0.0-1.0, default: 0.3)

#### Response

```json
{
  "query_smiles": "c1ccccc1",
  "num_matches": 5,
  "matches": [
    {
      "cluster_id": 12345,
      "similarity": 0.875,
      "centroid_smiles": "c1ccncc1",
      "num_members": 150,
      "example_smiles": [
        "c1ccncc1",
        "c1cnccc1",
        "c1cncc1"
      ]
    }
  ]
}
```

**Response Fields:**
- `query_smiles`: Original query ring SMILES
- `num_matches`: Number of matches found
- `matches`: Array of matched clusters
  - `cluster_id`: Unique ID of the cluster
  - `similarity`: Tanimoto similarity score (0.0-1.0)
  - `centroid_smiles`: Representative SMILES of the cluster
  - `num_members`: Number of rings in this cluster
  - `example_smiles`: Up to 5 example SMILES from the cluster

#### Example cURL Request

```bash
curl -X POST "http://localhost:8000/molecule/bioisostere/scan" \
  -H "Content-Type: application/json" \
  -d '{
    "ring_smiles": "c1ccc2ccccc2c1",
    "top_k": 10,
    "min_similarity": 0.4
  }'
```

### 2. Get Cluster Details

**GET** `/molecule/bioisostere/cluster/{cluster_id}`

Retrieve detailed information about a specific cluster.

#### Path Parameters

- `cluster_id` (integer): ID of the cluster

#### Response

```json
{
  "cluster_id": 12345,
  "num_members": 150,
  "member_smiles": [
    "c1ccncc1",
    "c1cnccc1",
    "c1cncc1",
    "..."
  ]
}
```

#### Example cURL Request

```bash
curl -X GET "http://localhost:8000/molecule/bioisostere/cluster/27678"
```

## Python Usage Examples

### Basic Scan

```python
import requests

# Scan for bio-isosteres of naphthalene
response = requests.post(
    "http://localhost:8000/molecule/bioisostere/scan",
    json={
        "ring_smiles": "c1ccc2ccccc2c1",
        "top_k": 10,
        "min_similarity": 0.4
    }
)

data = response.json()
print(f"Found {data['num_matches']} matches")

for i, match in enumerate(data['matches'], 1):
    print(f"{i}. Similarity: {match['similarity']:.3f}")
    print(f"   Representative: {match['centroid_smiles']}")
    print(f"   Cluster size: {match['num_members']} rings")
```

### Direct Service Usage (within backend)

```python
from src.services.mcp.ring_scan import get_scanner

# Initialize scanner
scanner = get_scanner()

# Scan for bio-isosteres
results = scanner.scan(
    query_smiles="c1ccncc1",  # pyridine
    top_k=20,
    min_similarity=0.3
)

for result in results:
    print(f"Cluster {result['cluster_id']}: {result['similarity']:.3f}")
    print(f"  Examples: {result['example_smiles'][:3]}")
```

### Get Full Cluster Information

```python
from src.services.mcp.ring_scan import get_scanner

scanner = get_scanner()

# Get all members of a cluster
cluster_info = scanner.get_cluster_info(27678)
print(f"Cluster contains {cluster_info['num_members']} rings")
print(f"All SMILES: {cluster_info['member_smiles']}")
```

## Test Results

The scanner has been tested with various ring structures:

| Query Ring | SMILES | Matches Found | Top Similarity |
|------------|--------|---------------|----------------|
| Benzene | `c1ccccc1` | 0 | N/A |
| Naphthalene | `c1ccc2ccccc2c1` | 10 | 0.875 |
| Pyridine | `c1ccncc1` | 2 | 0.353 |
| Furan | `c1ccoc1` | 0 | N/A |
| Imidazole | `c1[nH]cnc1` | 0 | N/A |
| Cyclohexane | `C1CCCCC1` | 2 | 0.333 |
| Piperidine | `C1CCNCC1` | 7 | 0.400 |

## Algorithm Details

### Fingerprint Generation

1. **Type**: Morgan fingerprints (circular fingerprints)
2. **Radius**: 2 (equivalent to ECFP4)
3. **Size**: 2048 bits
4. **Packing**: 8 bits per byte for efficient storage

### Similarity Calculation

The algorithm uses **Tanimoto coefficient** (Jaccard similarity) for binary fingerprints:

```
Tanimoto(A, B) = |A ∩ B| / |A ∪ B|
               = intersection / (cardinality_A + cardinality_B - intersection)
```

Where:
- `A` and `B` are binary fingerprints
- `|A ∩ B|` is the number of bits set in both fingerprints
- `|A ∪ B|` is the total number of bits set in either fingerprint

### Performance

- **Initialization**: ~2-3 seconds (loads ~34k centroids and ~4M SMILES)
- **Per-query scan**: ~1-2 seconds (compares against all centroids)
- **Memory usage**: ~500MB for loaded data structures

### Clustering Method

The database was pre-clustered using:
- **Algorithm**: BIRCH (Balanced Iterative Reducing and Clustering using Hierarchies)
- **Similarity metric**: Tanimoto similarity on Morgan fingerprints
- **Threshold**: 0.3 (minimum similarity for cluster membership)
- **Branching factor**: 254
- **Result**: 34,389 clusters from 3,931,782 rings

## Use Cases

### 1. Scaffold Hopping
Find alternative ring scaffolds with similar properties:
```python
# Find replacements for a benzimidazole core
results = scanner.scan("c1ccc2[nH]cnc2c1", top_k=20)
```

### 2. Lead Optimization
Explore bioisosteric replacements to improve ADME properties:
```python
# Find alternatives to pyridine to modify pKa
results = scanner.scan("c1ccncc1", min_similarity=0.4)
```

### 3. Patent Circumvention
Identify structurally distinct but similar rings:
```python
# Find rings similar to a patented scaffold
results = scanner.scan("c1ccc2c(c1)ccn3cccc23", top_k=30)
```

### 4. Library Design
Generate diverse compound libraries around a core structure:
```python
# Get diverse analogs of indole
results = scanner.scan("c1ccc2[nH]ccc2c1", min_similarity=0.5, top_k=50)
```

## Limitations

1. **Query must be a ring structure**: Non-cyclic structures are not supported
2. **Small rings may have few matches**: Simple 5-6 membered rings may return no matches above the similarity threshold
3. **Threshold sensitivity**: Default 0.3 threshold may miss slightly similar structures; adjust based on needs
4. **Centroid representation**: Returned examples are cluster members, not exact centroids

## File Structure

```
backend/src/
├── services/
│   └── mcp/
│       └── ring_scan/                 # Ring scan MCP tool
│           ├── __init__.py            # Module exports
│           ├── tools/
│           │   ├── __init__.py
│           │   └── bioisostere_scanner.py  # Main scanner service
│           └── db/
│               ├── rings.smi           # All ring SMILES (~4M rings)
│               ├── rings.npy           # Packed fingerprints
│               └── clusters/
│                   ├── cluster-centroids-packed.pkl  # 34k centroids
│                   ├── clusters.pkl    # Cluster assignments
│                   └── config.json     # Clustering configuration
├── controllers/
│   └── bioisostere_controller.py      # API controller
├── api/
│   ├── models/molecule.py             # Request/response models
│   └── endpoints/molecule.py          # API endpoints
└── test_bioisostere_scanner.py        # Test script
```

## Running Tests

```bash
cd backend
python test_bioisostere_scanner.py
```

## Dependencies

- `rdkit>=2025.9.3`: Cheminformatics toolkit for SMILES and fingerprints
- `numpy>=2.3.5`: Array operations and bit manipulation
- `bblean>=0.6.1b0`: Clustering utilities (optional)
- `fastapi>=0.104.0`: API framework
- `pydantic>=2.5.0`: Data validation

## References

1. Ertl, P. "Database of 4 million medicinal chemistry relevant ring systems" J. Chem. Inf. Model. 2024
2. Rogers, D. & Hahn, M. "Extended-Connectivity Fingerprints" J. Chem. Inf. Model. 2010
3. Zhang, T., Ramakrishnan, R., & Livny, M. "BIRCH: An efficient data clustering method for very large databases" SIGMOD 1996
