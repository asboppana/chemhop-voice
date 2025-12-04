"""
Pydantic models for MCP tool responses.

These models ensure type safety when the agent returns structured data from MCP tools.
The agent can call various MCP tools during a conversation, and we want to:
1. Type the responses properly (not just Dict[str, Any])
2. Separate natural language text from structured JSON data
3. Track which tools were called and what they returned

=============================================================================
RESPONSE TYPE PATTERNS
=============================================================================

The agent workflows produce 4 distinct response patterns:

1. MOLECULE ANNOTATION (Phase 1 MCP - annotate_molecule)
   - Input: User text (e.g., "analyze aspirin")
   - Output: SMILES + full annotation with functional groups, SVG, matches
   - Model: MoleculeAnnotation
   - ONE annotation per molecule
   
2. SMILES QUERY SET (Agent reasoning step)
   - Input: User text (e.g., "find bioisosteres for the benzene and pyridine rings")
   - Output: List of targets to query (agent extracts from annotation)
   - Model: SmilesQuerySet with multiple SmilesQueryTarget entries
   - Note: Each target will get its OWN MCP tool call

3. BIOISOSTERE RESULTS (Phase 1 MCP - scan_for_bioisosteres)
   - Input: Single SMILES (one ring/fragment)
   - Output: Many bioisostere candidates with similarity scores
   - Model: BioisostereScanResult
   - ONE result per target (if 3 targets requested, 3 separate scan calls, 3 results)
   
4. ADMET PROPERTIES (Phase 2 MCP - predict_admet_properties)
   - Input: Single SMILES
   - Output: Dictionary of scalar property predictions
   - Model: ADMETPropertyResult
   - ONE result per SMILES queried

=============================================================================
KEY PRINCIPLE: ONE MCP CALL PER TARGET
=============================================================================

When user requests bioisosteres/properties for multiple groups:
- Agent extracts SMILES for EACH group → SmilesQuerySet with N targets
- Agent makes N SEPARATE MCP tool calls (one per target)
- Agent emits N SEPARATE structured JSON results (one per call)

Example: "find bioisosteres for benzene and pyridine"
→ SmilesQuerySet with 2 targets
→ scan_for_bioisosteres("c1ccccc1") → BioisostereStructuredMessage #1
→ scan_for_bioisosteres("c1ccncc1") → BioisostereStructuredMessage #2

=============================================================================
MCP SERVERS
=============================================================================

Phase 1 MCP (drugdiscovery_mcp):
- annotate_molecule: SMILES → MoleculeAnnotation
- get_smiles_from_name: chemical_name → SmilesFromName
- convert_identifier_to_smiles: identifier → IdentifierToSmiles
- check_molecule_patent: SMILES → PatentCheckResult
- scan_for_bioisosteres: SMILES → BioisostereScanResult (CALL ONCE PER TARGET)

Phase 2 MCP (admet_mcp):
- predict_admet_properties: SMILES → ADMETPropertyResult (CALL ONCE PER SMILES)

=============================================================================
"""
from typing import List, Optional, Dict, Any, Union
from pydantic import BaseModel, Field


# ============================================================================
# Annotate Molecule Models
# ============================================================================

class TrivialName(BaseModel):
    """Functional group pattern information."""
    name: str = Field(..., description="Common name of the functional group")
    smarts: str = Field(..., description="SMARTS pattern for the functional group")
    group: str = Field(..., description="Category (e.g., 'cyclic', 'acyclic', 'functional_group')")
    bonds: int = Field(..., description="Number of bonds in the pattern")
    hierarchy: Optional[str] = Field(None, description="Hierarchical relationship to other patterns")


class AnnotationMatch(BaseModel):
    """A single functional group match in the molecule."""
    atom_indices: List[int] = Field(..., description="Indices of atoms matching the pattern")
    trivial_name: TrivialName = Field(..., description="Pattern information")


class MoleculeAnnotation(BaseModel):
    """Complete annotation of a molecule with functional groups."""
    name: Optional[str] = Field(None, description="Molecule name (if available)")
    svg: str = Field(..., description="SVG image representation of the molecule")
    smiles: str = Field(..., description="Canonical SMILES string")
    matches: List[AnnotationMatch] = Field(default_factory=list, description="List of matched functional groups")
    error: Optional[str] = Field(None, description="Error message if annotation failed")


# ============================================================================
# Identifier Conversion Models
# ============================================================================

class SmilesFromName(BaseModel):
    """SMILES string retrieved from chemical name."""
    smiles: Optional[str] = Field(None, description="Canonical SMILES string")
    error: Optional[str] = Field(None, description="Error message if lookup failed")


class IdentifierToSmiles(BaseModel):
    """Conversion of chemical identifier to SMILES."""
    identifier: str = Field(..., description="Original input identifier")
    smiles_list: List[str] = Field(default_factory=list, description="List of SMILES strings")
    count: int = Field(0, description="Number of SMILES strings returned")
    error: Optional[str] = Field(None, description="Error message if conversion failed")


# ============================================================================
# Patent Search Models
# ============================================================================

class PatentInfo(BaseModel):
    """Individual patent information."""
    patent_id: str = Field(..., description="Patent ID")
    title: Optional[str] = Field(None, description="Patent title")
    publication_date: Optional[str] = Field(None, description="Publication date")
    url: Optional[str] = Field(None, description="URL to patent document")


class PatentCheckResult(BaseModel):
    """Patent status check result."""
    smiles: str = Field(..., description="SMILES string of the molecule")
    is_patented: Optional[bool] = Field(None, description="Whether molecule appears in patents")
    fto_status: Optional[str] = Field(None, description="Freedom to operate status")
    confidence: Optional[str] = Field(None, description="Confidence level of the result")
    patents: List[PatentInfo] = Field(default_factory=list, description="List of relevant patents")
    error: Optional[str] = Field(None, description="Error message if check failed")


# ============================================================================
# Bioisostere Scanning Models
# ============================================================================

class MolecularDescriptors(BaseModel):
    """Physicochemical properties of a molecule."""
    logp: Optional[float] = Field(None, description="Lipophilicity (logP)")
    tpsa: Optional[float] = Field(None, description="Topological polar surface area (Ų)")
    h_donors: Optional[int] = Field(None, description="Number of hydrogen bond donors")
    h_acceptors: Optional[int] = Field(None, description="Number of hydrogen bond acceptors")
    num_rings: Optional[int] = Field(None, description="Number of rings")
    num_aromatic_rings: Optional[int] = Field(None, description="Number of aromatic rings")
    molecular_weight: Optional[float] = Field(None, description="Molecular weight (g/mol)")


class DeltaProperties(BaseModel):
    """Property differences from query molecule."""
    delta_logp: Optional[float] = Field(None, description="Difference in logP")
    delta_tpsa: Optional[float] = Field(None, description="Difference in TPSA")
    delta_h_donors: Optional[int] = Field(None, description="Difference in H-donors")
    delta_h_acceptors: Optional[int] = Field(None, description="Difference in H-acceptors")


class BioisostereResult(BaseModel):
    """A single bio-isosteric replacement candidate."""
    source: str = Field(..., description="Data source (ertl, chemspace, or clusters)")
    centroid_smiles: str = Field(..., description="SMILES of the bio-isostere")
    similarity: float = Field(..., ge=0.0, le=1.0, description="Structural similarity score")
    bio_isostere_score: float = Field(..., ge=0.0, le=1.0, description="Combined bio-isostere score")
    pharmacophore_similarity: Optional[float] = Field(None, ge=0.0, le=1.0, description="2D pharmacophore similarity")
    topology_similarity: Optional[float] = Field(None, ge=0.0, le=1.0, description="Ring topology similarity")
    descriptors: MolecularDescriptors = Field(..., description="Physicochemical properties")
    delta_properties: Optional[DeltaProperties] = Field(None, description="Property differences from query")
    cluster_id: Optional[str] = Field(None, description="ID of the cluster (for cluster results)")
    num_members: Optional[int] = Field(None, description="Number of members in cluster")
    example_smiles: Optional[List[str]] = Field(None, description="Example SMILES from cluster")


class BioisostereScanResult(BaseModel):
    """
    Results from bio-isostere scanning for a SINGLE query.
    
    Each scan is for one SMILES query. If user asked for multiple groups,
    there will be multiple BioisostereScanResult objects - one per group.
    """
    query_smiles: str = Field(..., description="Input query SMILES")
    source_group: Optional[str] = Field(None, description="Name of the group this came from (e.g., 'benzene ring')")
    atom_indices: Optional[List[int]] = Field(None, description="Atom indices from original molecule")
    num_results: int = Field(0, description="Number of bio-isosteres found")
    results: List[BioisostereResult] = Field(default_factory=list, description="List of bio-isosteric replacements")
    error: Optional[str] = Field(None, description="Error message if scan failed")


# ============================================================================
# ADMET Property Prediction Models (Phase 2 MCP - admet_mcp)
# ============================================================================

class ADMETPropertyResult(BaseModel):
    """
    ADMET property predictions for a molecule.
    
    The predictions dict contains ~50+ properties including:
    - Absorption: Caco2_Wang, Solubility_AqSolDB, HIA_Hou, Pgp_Broccatelli, etc.
    - Distribution: BBB_Martins, PPB_AZ, VDss_Lombardo, etc.
    - Metabolism: CYP1A2/2C9/2C19/2D6/3A4 inhibition, etc.
    - Excretion: Half_Life_Obach, Clearance_Hepatocyte, etc.
    - Toxicity: hERG, AMES, DILI, Skin_Reaction, etc.
    """
    smiles: str = Field(..., description="Input SMILES string")
    predictions: Dict[str, float] = Field(default_factory=dict, description="Property name → predicted value")
    error: Optional[str] = Field(None, description="Error message if prediction failed")


# ============================================================================
# SMILES Query Set (Agent reasoning output - before tool calls)
# ============================================================================

class SmilesQueryTarget(BaseModel):
    """
    A single SMILES target extracted from user intent.
    
    Each target represents one group/fragment that the user wants to query.
    The agent will make a SEPARATE MCP tool call for each target.
    """
    smiles: str = Field(..., description="SMILES string to query")
    source_group: str = Field(..., description="Name of the annotated group (e.g., 'benzene ring', 'pyridine ring')")
    atom_indices: List[int] = Field(default_factory=list, description="Atom indices from the original molecule annotation")


class SmilesQuerySet(BaseModel):
    """
    Set of SMILES strings to query, extracted from user intent.
    
    When user says "find bioisosteres for the benzene and pyridine rings",
    the agent:
    1. Extracts SMILES for each referenced group from the annotation
    2. Creates a target entry for EACH group
    3. Will make SEPARATE MCP tool calls for each target
    
    There is NO boolean logic - each target is queried independently.
    """
    original_molecule_smiles: Optional[str] = Field(None, description="Parent molecule SMILES (if applicable)")
    targets: List[SmilesQueryTarget] = Field(default_factory=list, description="List of targets to query - one MCP call per target")
    query_type: str = Field(..., description="Type of query: 'bioisostere', 'admet', or 'combined'")
    reasoning: Optional[str] = Field(None, description="Agent's reasoning for which groups were selected")


# ============================================================================
# Combined Multi-Query Results
# ============================================================================

class BioisostereQueryResult(BaseModel):
    """Bioisostere results for a single query in a multi-query operation."""
    query: SmilesQueryTarget = Field(..., description="The query that was made")
    scan_result: BioisostereScanResult = Field(..., description="Results from the bioisostere scan")


class ADMETQueryResult(BaseModel):
    """ADMET results for a single query in a multi-query operation."""
    query: SmilesQueryTarget = Field(..., description="The query that was made")
    admet_result: ADMETPropertyResult = Field(..., description="ADMET predictions")


class MultiQueryBioisostereResult(BaseModel):
    """
    Results from querying bioisosteres for multiple targets.
    
    When user says "find bioisosteres for the benzene and pyridine rings",
    this contains results for each individual target queried.
    Each target gets its own MCP call and its own result entry.
    """
    query_set: SmilesQuerySet = Field(..., description="The original query set")
    results: List[BioisostereQueryResult] = Field(default_factory=list, description="One result per target queried")
    total_bioisosteres_found: int = Field(0, description="Total across all queries")


class MultiQueryADMETResult(BaseModel):
    """
    Results from predicting ADMET properties for multiple SMILES.
    """
    query_set: SmilesQuerySet = Field(..., description="The original query set")
    results: List[ADMETQueryResult] = Field(default_factory=list, description="Results per query")


# ============================================================================
# Agent Response Models
# ============================================================================

class MCPToolCall(BaseModel):
    """Represents a single MCP tool invocation and its result."""
    tool_name: str = Field(..., description="Name of the MCP tool called")
    server: str = Field(default="phase1", description="Which MCP server: 'phase1' or 'phase2'")
    tool_args: Dict[str, Any] = Field(..., description="Arguments passed to the tool")
    result: Union[
        # Phase 1 MCP tools
        MoleculeAnnotation,
        SmilesFromName,
        IdentifierToSmiles,
        PatentCheckResult,
        BioisostereScanResult,
        # Phase 2 MCP tools
        ADMETPropertyResult,
        # Fallback
        Dict[str, Any]
    ] = Field(..., description="Result from the tool call")


class AgentChatResponse(BaseModel):
    """
    Complete response from the chat agent.
    
    The response contains:
    - response: Natural language text (cleaned of any embedded JSON)
    - status: "success" or "error"
    - tool_calls: List of MCP tools that were called and their results
    - structured_data: Any embedded structured data extracted from the response
    - response_type: Indicates what type of structured data is present
    """
    response: str = Field(..., description="Natural language response from the agent")
    status: str = Field(default="success", description="Status of the request")
    response_type: Optional[str] = Field(
        None, 
        description="Type of response: 'annotation', 'query_set', 'bioisosteres', 'admet', 'multi_bioisostere', 'multi_admet', 'text_only'"
    )
    tool_calls: Optional[List[MCPToolCall]] = Field(None, description="MCP tools called during processing")
    structured_data: Optional[Dict[str, Any]] = Field(None, description="Any structured data extracted from the response")
    
    # Typed structured outputs (mutually exclusive based on response_type)
    molecule_data: Optional[MoleculeAnnotation] = Field(None, description="Present when response_type='annotation'")
    query_set_data: Optional[SmilesQuerySet] = Field(None, description="Present when response_type='query_set'")
    bioisostere_data: Optional[BioisostereScanResult] = Field(None, description="Present when response_type='bioisosteres'")
    admet_data: Optional[ADMETPropertyResult] = Field(None, description="Present when response_type='admet'")
    multi_bioisostere_data: Optional[MultiQueryBioisostereResult] = Field(None, description="Present when response_type='multi_bioisostere'")
    multi_admet_data: Optional[MultiQueryADMETResult] = Field(None, description="Present when response_type='multi_admet'")
    
    class Config:
        json_schema_extra = {
            "example": {
                "response": "I've annotated aspirin and found 2 functional groups...",
                "status": "success",
                "response_type": "annotation",
                "tool_calls": [
                    {
                        "tool_name": "annotate_molecule",
                        "server": "phase1",
                        "tool_args": {"smiles": "CC(=O)Oc1ccccc1C(=O)O"},
                        "result": {
                            "name": "aspirin",
                            "svg": "<svg>...</svg>",
                            "smiles": "CC(=O)Oc1ccccc1C(=O)O",
                            "matches": []
                        }
                    }
                ],
                "molecule_data": {
                    "name": "aspirin",
                    "svg": "<svg>...</svg>",
                    "smiles": "CC(=O)Oc1ccccc1C(=O)O",
                    "matches": []
                }
            }
        }


# ============================================================================
# Structured Message Types (for embedded JSON in response text)
# ============================================================================

class MoleculeStructuredMessage(BaseModel):
    """
    Structured message for molecule annotation data.
    Emitted when agent annotates a molecule.
    """
    type: str = Field(default="molecule_structured", description="Message type identifier")
    smiles: str = Field(..., description="SMILES string")
    annotation: MoleculeAnnotation = Field(..., description="Full annotation object")


class QuerySetStructuredMessage(BaseModel):
    """
    Structured message for SMILES query sets.
    Emitted when agent extracts SMILES to query from user intent.
    """
    type: str = Field(default="query_set_structured", description="Message type identifier")
    query_set: SmilesQuerySet = Field(..., description="Set of SMILES to query")


class BioisostereStructuredMessage(BaseModel):
    """
    Structured message for bioisostere results.
    
    Emitted ONCE for EACH target that was scanned.
    If user asked for bioisosteres for 3 groups, there will be 3 of these messages.
    """
    type: str = Field(default="bioisostere_structured", description="Message type identifier")
    query_smiles: str = Field(..., description="SMILES that was queried")
    source_group: str = Field(..., description="Name of the group this came from (e.g., 'benzene ring')")
    atom_indices: List[int] = Field(default_factory=list, description="Atom indices from the original molecule")
    results: List[BioisostereResult] = Field(..., description="List of bioisostere candidates found")


class ADMETStructuredMessage(BaseModel):
    """
    Structured message for ADMET property predictions.
    
    Emitted ONCE for EACH molecule/fragment that was predicted.
    """
    type: str = Field(default="admet_structured", description="Message type identifier")
    smiles: str = Field(..., description="SMILES that was predicted")
    source_group: Optional[str] = Field(None, description="Name of the group if this is a fragment (e.g., 'benzene ring')")
    predictions: Dict[str, float] = Field(..., description="Property name → predicted value")


# Union type for all structured messages
StructuredMessage = Union[
    MoleculeStructuredMessage,
    QuerySetStructuredMessage,
    BioisostereStructuredMessage,
    ADMETStructuredMessage
]