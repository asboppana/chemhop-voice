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
