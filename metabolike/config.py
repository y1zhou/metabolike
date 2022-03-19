import logging
from pathlib import Path
from typing import Optional

import yaml
from pydantic import BaseModel

logger = logging.getLogger(__name__)


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


class BrendaConfig(BaseModel):
    brenda_file: str
    bkms_react_file: str


class Config(BaseModel):
    database: DatabaseConfig
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
