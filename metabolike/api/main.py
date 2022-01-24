import logging

from fastapi import FastAPI

logger = logging.getLogger(__name__)
app = FastAPI()


@app.get("/")
def read_root():
    return {"Hello": "World"}


# @app.get("/super_pathways")
# def read_super_pathways():
#     return super_pathways
