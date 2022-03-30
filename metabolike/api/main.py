import logging

from fastapi import FastAPI

from .core.db import db
from .routers import pathways

logger = logging.getLogger(__name__)
app = FastAPI()

app.include_router(pathways.router)


@app.get("/")
def read_root():
    return {"Hello": "World"}


# Shutdown event handler
@app.on_event("shutdown")
def shutdown():
    db.close()
