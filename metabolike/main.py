#!/usr/bin/env python
import logging
from pathlib import Path
from typing import Optional

import typer
import uvicorn
import yaml

from metabolike.config import Config
from metabolike.db.metacyc import MetaDB
from metabolike.parser.brenda import parse_brenda
from metabolike.parser.metacyc import Metacyc

logging.basicConfig(
    format="[%(levelname)s] %(asctime)s - %(name)s:%(lineno)s:%(funcName)s - %(message)s",
    level=logging.INFO,
)
logger = logging.getLogger(__name__)

app = typer.Typer()


def load_config(config_file: str) -> Config:
    logger.info(f"Parsing config file {config_file}")
    conf_file = Path(config_file).expanduser().resolve()
    with conf_file.open("r") as f:
        conf_data = yaml.safe_load(f)

    # Validate config file using pydantic schema
    conf = Config(**conf_data)
    return conf


@app.command()
def setup(
    config_file: str = typer.Option(
        ..., "--config", "-c", help="Path to the configuration file."
    ),
    database: Optional[str] = typer.Option(
        None,
        help="Name of the database. Use MetaID from SBML file if not given.",
    ),
    create_db: bool = typer.Option(
        True, help="When database creation is not allowed, set this to False."
    ),
    drop_if_exists: bool = typer.Option(
        False, "--drop-if-exists", "-f", help="Drop the database if it already exists."
    ),
):
    conf = load_config(config_file)

    logger.info("Connecting to neo4j database")
    db = MetaDB(**conf.database.dict())
    meta = Metacyc(db, **conf.metacyc.dict(), db_name=database)

    logger.info("Setting up database using MetaCyc data")
    meta.setup(create_db=create_db, force=drop_if_exists)

    if conf.brenda:
        logger.info("Reading BRENDA text file")
        brenda = parse_brenda(conf.brenda.brenda_file)

    meta.db.close()


@app.command()
def serve():
    """Serves an API to the graph database.

    The command uses the following environment variables to configure the
    database connection: METABOLIKE_DB_URI, METABOLIKE_DB_USER,
    METABOLIKE_DB_PASSWORD, and METABOLIKE_DB_NAME.

    Alternatively, you may save the variables in a file named "metabolike.env"
    in the current directory.
    """
    uvicorn.run(
        "metabolike.api.main:app",
        host="0.0.0.0",
        port=8000,
        log_level="info",
        proxy_headers=True,
    )


if __name__ == "__main__":
    app()
