# SureChembl Tool for patent search

from typing import Dict, Any
from rdkit import Chem
import requests
import logging

logger = logging.getLogger(__name__)

class SureChEMBLPatentTool:

    def __init__(self):
        self.chembl_base = "https://www.ebi.ac.uk/chembl/api/data"
        self.surechembl_base = "https://www.surechembl.org/chemical"

    def check_patent(self, smiles: str) -> Dict[str, Any]:
        """Check is molecule appears in patent literature returns FTO status and patent details. """

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {smiles}")
        
        canonical_smiles = Chem.MolToSmiles(mol)
        inchi_key = Chem.MolToInchiKey(mol)
        
        result = self._search_chembl(canonical_smiles, inchi_key)
        
        return {
            "smiles": canonical_smiles,
            "inchi_key": inchi_key,
            "is_patented": result["found"],
            "patent_count": result["patent_count"],
            "surechembl_id": result.get("surechembl_id"),
            "patents": result.get("patents", []),
            "surechembl_url": f"{self.surechembl_base}/{inchi_key}" if result["found"] else None,
            "fto_status": "POTENTIAL_CONFLICT" if result["found"] else "LIKELY_CLEAR",
            "confidence": "high" if result["found"] else "medium",
            "note": "Not found in patent database." if not result["found"] else "Found in patent literature."
        }
    
    def _search_chembl(self, smiles: str, inchi_key: str) -> Dict[str, Any]:
        """Search ChEMBL/SureChEMBL for compound."""
        result = {
            "found": False,
            "patent_count": 0,
            "patents": [],
            "surechembl_id": None
        }
        
        try:
            url = f"{self.chembl_base}/molecule/{inchi_key}.json"
            response = requests.get(url, timeout=10)
            
            if response.status_code == 200:
                data = response.json()
                
                if "molecule_properties" in data:
                    result["found"] = True
                    result["surechembl_id"] = data.get("molecule_chembl_id")
                    self._extract_patent_info(data, result)
                    
                logger.info(f"Found {inchi_key} in ChEMBL/SureChEMBL")
            else:
                # result = self._search_by_similarity(smiles)
                logger.info(f"No exact match for {inchi_key}")
                
        except requests.exceptions.RequestException as e:
            logger.warning(f"ChEMBL search failed: {e}")
        except Exception as e:
            logger.error(f"Error in ChEMBL search: {e}")
        
        return result



#    def _search_by_similarity(self, smiles: str) -> Dict[str, Any]:
#         """Fallback: search by similarity if exact match fails."""
#         result = {"found": False, "patent_count": 0, "patents": [], "surechembl_id": None}
        
#         try:
#             url = f"{self.chembl_base}/similarity/{smiles}/70.json"
#             response = requests.get(url, timeout=15)
            
#             if response.status_code == 200:
#                 data = response.json()
#                 molecules = data.get("molecules", [])
                
#                 if molecules:
#                     for mol in molecules[:5]:
#                         if mol.get("molecule_type") == "Small molecule":
#                             result["found"] = True
#                             result["patent_count"] += 1
                    
#                     logger.info(f"Found {len(molecules)} similar molecules")
#         except Exception as e:
#             logger.error(f"Similarity search failed: {e}")
        
#         return result
    def _extract_patent_info(self, data: Dict, result: Dict) -> None:
        """Extract patent info from ChEMBL response."""
        if "cross_references" in data:
            for ref in data.get("cross_references", []):
                if ref.get("xref_src") in ["PubChem", "SureChEMBL"]:
                    result["patent_count"] += 1
                    result["patents"].append({
                        "source": ref.get("xref_src"),
                        "id": ref.get("xref_id")
                    })
    