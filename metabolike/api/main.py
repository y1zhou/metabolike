import atexit
import logging

from fastapi import FastAPI
from metabolike.api.config import db_conf
from metabolike.db.metacyc import MetaDB

logger = logging.getLogger(__name__)
app = FastAPI()

db = MetaDB(
    db_conf.metabolike_db_uri,
    db_conf.metabolike_db_user,
    db_conf.metabolike_db_password,
)
db.use_database(db_conf.metabolike_db_name)


@app.get("/")
def read_root():
    return {"Hello": "World"}


@app.get("/pathway/{pathway_id}")
def get_pathway(pathway_id: str):
    nodes, edges = db.get_view_of_pathway(pathway_id)
    return {"nodes": [x[0] for x in nodes], "edges": [x[0] for x in edges]}


def exit_application():
    db.close()


atexit.register(exit_application)
