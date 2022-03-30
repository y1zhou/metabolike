import logging

from fastapi import APIRouter

from ..core.db import db
from ..schemas.pathways import Pathways

logger = logging.getLogger(__name__)

router = APIRouter()


@router.get("/pathway/{pathway_id}", tags=["pathways"])
def get_pathway(pathway_id: str):
    nodes, edges = db.get_view_of_pathway(pathway_id)
    return {"nodes": nodes, "edges": edges}


@router.post("/pathways/reactions/", tags=["pathways"])
async def get_pathway_reactions(pw: Pathways):
    logger.debug(f"Request for reactions in pathways: {pw.pathway_ids}")
    df, rxn_genes = db.get_fba_info_of_pathways(pw.pathway_ids)
    return {"reactions": df.to_dict(), "genes": rxn_genes}
