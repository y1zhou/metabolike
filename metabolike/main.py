#!/usr/bin/env python
import logging

import typer
import uvicorn

from .config import load_config
from .db import MetacycClient
from .parser import MetacycParser
from .parser.brenda import parse_brenda

logging.basicConfig(
    format="[%(levelname)s] %(asctime)s - %(name)s:%(lineno)s:%(funcName)s - %(message)s",
    level=logging.INFO,
)
logger = logging.getLogger(__name__)

app = typer.Typer()


@app.command()
def setup(
    config_file: str = typer.Option(
        ..., "--config", "-c", help="Path to the configuration file."
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
    db = MetacycClient(
        **conf.neo4j.dict(include={"uri", "database"}),
        neo4j_user=conf.neo4j.user,
        neo4j_password=conf.neo4j.password.get_secret_value(),
        create_db=create_db,
        drop_if_exists=drop_if_exists,
    )
    meta = MetacycParser(
        db,
        **conf.metacyc.dict(),
    )

    logger.info("Setting up database using MetaCyc data")
    meta.setup()

    # if conf.brenda:
    #     logger.info("Reading BRENDA text file")
    #     all_ecs = db.get_all_ec_numbers()
    # brenda = parse_brenda(conf.brenda.brenda_file, ec_nums=all_ecs, cache=use_cache)

    db.close()


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
