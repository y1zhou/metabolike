"""Schema for the YAML config files passed to `metabolike setup`.

The main sections are [config.Neo4jConfig][] and [config.MetacycConfig][].
"""
import logging
from pathlib import Path
from typing import Optional

import yaml
from pydantic import BaseModel, SecretStr

logger = logging.getLogger(__name__)


class Neo4jConfig(BaseModel):
    """Config schema for connecting to a Neo4j instance.

    See [db.neo4j.Neo4jClient][] for details. Valid keys are `uri`, `user`, `password`, and
    `database`.
    """

    uri: str
    user: str
    password: SecretStr
    database: str = "neo4j"


class MetacycConfig(BaseModel):
    """Config schema for parsing MetaCyc data files into a graph.

    Valid keys include `sbml`, `reactions`, `atom_mapping`, `pathways`, `compounds`,
    `publications`, and `classes`. All values should be paths to local files. Only the `sbml` key
    is required.
    """

    sbml: str
    reactions: Optional[str]
    atom_mapping: Optional[str]
    pathways: Optional[str]
    compounds: Optional[str]
    publications: Optional[str]
    classes: Optional[str]


class BrendaConfig(BaseModel):
    brenda_file: str
    bkms_react_file: str


class Config(BaseModel):
    neo4j: Neo4jConfig
    metacyc: MetacycConfig
    brenda: Optional[BrendaConfig]


def load_config(config_file: str) -> Config:
    logger.info(f"Parsing config file {config_file}")
    conf_file = Path(config_file).expanduser().resolve()
    with conf_file.open("r") as f:
        conf_data = yaml.safe_load(f)

    # Validate config file using pydantic schema
    conf = Config(**conf_data)
    return conf
