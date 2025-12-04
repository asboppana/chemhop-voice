import axios from 'axios';

const API_BASE_URL = import.meta.env.VITE_API_BASE_URL || 'http://localhost:8000';

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

export interface BioisostereMatch {
  source: string;
  similarity: number;
  centroid_smiles: string;
  bio_isostere_score: number;
  pharmacophore_similarity: number;
  topology_similarity: number;
  cluster_id?: number;
  num_members: number;
  example_smiles: string[];
  descriptors: Record<string, any>;
  delta_properties: Record<string, number>;
}

export interface BioisostereScanResponse {
  query_smiles: string;
  num_matches: number;
  matches: BioisostereMatch[];
}

export const scanBioisosteres = async (
  ringSmiles: string,
  topK: number = 20,
  minSimilarity: number = 0.3
): Promise<BioisostereScanResponse> => {
  const response = await axios.post(`${API_BASE_URL}/api/v1/molecule/bioisostere/scan`, {
    ring_smiles: ringSmiles,
    top_k: topK,
    min_similarity: minSimilarity
  });
  return response.data;
};

export interface ExtractSubstructureResponse {
  substructure_smiles: string;
  parent_smiles: string;
}

export const extractSubstructure = async (
  smiles: string,
  atomIndices: number[]
): Promise<ExtractSubstructureResponse> => {
  const response = await axios.post(`${API_BASE_URL}/api/v1/molecule/extract-substructure`, {
    smiles,
    atom_indices: atomIndices
  });
  return response.data;
};

export interface GenerateSvgResponse {
  svg: string;
  smiles: string;
}

export const generateSvg = async (
  smiles: string,
  width: number = 200,
  height: number = 200
): Promise<GenerateSvgResponse> => {
  const response = await axios.post(`${API_BASE_URL}/api/v1/molecule/generate-svg`, {
    smiles,
    width,
    height
  });
  return response.data;
};

export interface ReplaceSubstructureResponse {
  original_smiles: string;
  replacement_smiles_list: string[];
  num_generated: number;
}

export const replaceSubstructure = async (
  parentSmiles: string,
  sourcePatternSmiles: string,
  replacementSmiles: string,
  atomIndices: number[]
): Promise<ReplaceSubstructureResponse> => {
  const response = await axios.post(`${API_BASE_URL}/api/v1/molecule/replace-substructure`, {
    parent_smiles: parentSmiles,
    source_pattern_smiles: sourcePatternSmiles,
    replacement_smiles: replacementSmiles,
    atom_indices: atomIndices
  });
  return response.data;
};
