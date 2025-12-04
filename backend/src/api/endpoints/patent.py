"""Patent search endpoint."""
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from src.services.mcp.patent_search.surechembl_tool import SureChEMBLPatentTool

router = APIRouter()

class PatentCheckRequest(BaseModel):
    smiles: str

@router.post("/patent/check")
async def check_patent(request: PatentCheckRequest):
    """Quick patent check using SureChEMBL."""
    try:
        tool = SureChEMBLPatentTool()
        return tool.check_patent(request.smiles)
    except ValueError as e:
        raise HTTPException(status_code=422, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))