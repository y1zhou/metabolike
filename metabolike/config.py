from typing import Optional

from pydantic import BaseModel


class DatabaseConfig(BaseModel):
    uri: str
    neo4j_user: str
    neo4j_password: str


class MetacycConfig(BaseModel):
    sbml: str
    reactions: Optional[str]
    atom_mapping: Optional[str]
    pathways: Optional[str]
    compounds: Optional[str]
    publications: Optional[str]
    classes: Optional[str]


class Config(BaseModel):
    database: DatabaseConfig
    metacyc: MetacycConfig
