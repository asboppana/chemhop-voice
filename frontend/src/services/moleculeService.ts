import axios from 'axios';

const API_BASE_URL = import.meta.env.VITE_API_BASE_URL || 'http://localhost:8000';

export interface MoleculeAnalysisResponse {
  name: string;
  svg: string;
  matches: Array<{
    atom_indices: number[];
    trivial_name: {
      name: string;
      smarts: string;
      group: string;
      bonds: number;
      hierarchy: string;
      index: number;
    };
  }>;
  smiles: string;
}

export const analyzeMolecule = async (smiles: string): Promise<MoleculeAnalysisResponse> => {
  const response = await axios.post(`${API_BASE_URL}/api/v1/molecule/analyze`, {
    smiles
  });
  return response.data;
};
