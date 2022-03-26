import logging

from fastapi import FastAPI

from .core.db import db

logger = logging.getLogger(__name__)
app = FastAPI()


@app.get("/")
def read_root():
    return {"Hello": "World"}


@app.get("/pathway/{pathway_id}")
def get_pathway(pathway_id: str):
    nodes, edges = db.get_view_of_pathway(pathway_id)
    return {"nodes": nodes, "edges": edges}


@app.post("/pathways/reactions/")
async def get_pathway_reactions(pw: Pathways):
    logger.debug(f"Request for reactions in pathways: {pw.pathway_ids}")
    df, rxn_genes = db.get_fba_info_of_pathways(pw.pathway_ids)
    return {"reactions": df.to_dict(), "genes": rxn_genes}


# Shutdown event handler
@app.on_event("shutdown")
def shutdown():
    db.close()
