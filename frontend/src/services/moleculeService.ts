import axios from 'axios';

const API_BASE_URL = import.meta.env.VITE_API_BASE_URL || 'http://localhost:8000';

// ============================================================================
// Molecule Analysis Types (existing)
// ============================================================================

export interface MoleculeAnalysisResponse {
  name: string;
  svg: string;
  smiles: string;
  matches: Array<{
    atom_indices: number[];
    svg: string | null;
    trivial_name: {
      name: string;
      smarts: string;
      group: string;
      bonds: number;
      hierarchy: string;
      index: number;
    };
  }>;
}

export interface HighlightResponse {
  svg: string;
  smiles: string;
}

// ============================================================================
// Chat/Agent Response Types
// ============================================================================

export interface BioisostereResult {
  source: string;
  centroid_smiles: string;
  similarity: number;
  bio_isostere_score: number;
  pharmacophore_similarity?: number;
  topology_similarity?: number;
  descriptors: {
    logp?: number;
    tpsa?: number;
    h_donors?: number;
    h_acceptors?: number;
    num_rings?: number;
    num_aromatic_rings?: number;
    molecular_weight?: number;
  };
  delta_properties?: {
    delta_logp?: number;
    delta_tpsa?: number;
    delta_h_donors?: number;
    delta_h_acceptors?: number;
  };
}

export interface SmilesQueryTarget {
  smiles: string;
  source_group: string;
  atom_indices: number[];
}

export interface SmilesQuerySet {
  original_molecule_smiles?: string;
  targets: SmilesQueryTarget[];
  query_type: 'bioisostere' | 'admet' | 'combined';
  reasoning?: string;
}

export interface BioisostereScanResult {
  query_smiles: string;
  source_group?: string;
  atom_indices?: number[];
  num_results: number;
  results: BioisostereResult[];
  error?: string;
}

export interface ADMETPropertyResult {
  smiles: string;
  predictions: Record<string, number>;
  error?: string;
}

export type ChatResponseType = 
  | 'annotation' 
  | 'query_set' 
  | 'bioisosteres' 
  | 'admet' 
  | 'multi_bioisostere' 
  | 'multi_admet' 
  | 'text_only';

export interface AgentChatResponse {
  response: string;
  status: string;
  response_type?: ChatResponseType;
  tool_calls?: Array<{
    tool_name: string;
    server: string;
    tool_args: Record<string, any>;
    result: any;
  }>;
  structured_data?: {
    embedded_objects: any[];
  };
  // Typed structured outputs
  molecule_data?: MoleculeAnalysisResponse;
  query_set_data?: SmilesQuerySet;
  bioisostere_data?: BioisostereScanResult;
  admet_data?: ADMETPropertyResult;
  multi_bioisostere_data?: any;
  multi_admet_data?: any;
}

// ============================================================================
// API Functions
// ============================================================================

export const analyzeMolecule = async (smiles: string): Promise<MoleculeAnalysisResponse> => {
  const response = await axios.post(`${API_BASE_URL}/api/v1/molecule/analyze`, {
    smiles
  });
  return response.data;
};

export const highlightMolecule = async (smiles: string, atomIndices: number[]): Promise<HighlightResponse> => {
  const response = await axios.post(`${API_BASE_URL}/api/v1/molecule/highlight`, {
    smiles,
    atom_indices: atomIndices
  });
  return response.data;
};

export const chatWithAgent = async (
  messages: Array<{ role: 'user' | 'assistant'; content: string }>,
  model?: string
): Promise<AgentChatResponse> => {
  const response = await axios.post(`${API_BASE_URL}/chat`, {
    messages,
    model
  });
  return response.data;
};
