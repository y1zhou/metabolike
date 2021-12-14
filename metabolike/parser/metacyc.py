import logging
from pathlib import Path
from typing import Union

import libsbml
from metabolike import db

logger = logging.getLogger(__name__)


class Metacyc:
    def __init__(self, filepath: Union[str, Path], neo4j_driver: db.Neo4jDriver):
        self.neo4j_driver = neo4j_driver

        # Read SBML file
        self.filepath = Path(filepath).expanduser().resolve()
        self.doc = self.read_sbml()
        self.model: libsbml.Model = self.doc.getModel()
        self.db_name: str = self.model.getId().lower()

        # Setup Neo4j database
        self.setup_neo4j()

    def read_sbml(self) -> libsbml.SBMLDocument:
        reader = libsbml.SBMLReader()
        metacyc = reader.readSBML(self.filepath)

        for i in range(metacyc.getNumErrors()):
            err = metacyc.getError(i)
            logger.warning(
                "SBML reader raised %s at line %d, column %d: %s",
                err.getSeverityAsString(),
                err.getLine(),
                err.getColumn(),
                err.getMessage().replace("\n", " "),
            )

        return metacyc

    def setup_neo4j(self, **kwargs):
        # Create database
        try:
            db.create(self.neo4j_driver, self.db_name, **kwargs)
        except Exception as e:
            logger.fatal(f"Could not create database: {e}")
            raise

        # Set constraints
        with self.neo4j_driver.session(database=self.db_name) as session:
            logger.debug("Creating constraint for RDF nodes")
            r = session.run(
                """CREATE CONSTRAINT IF NOT EXISTS
                   ON (r:RDF) ASSERT r.uri IS UNIQUE;"""
            ).data()
            if r:
                logger.warning(f"Could not create constraint for RDF nodes: {r}")

            # Constraints automatically create indexes, so we don't need to
            # create them manually.
            for label in ["Pathway", "Compartment", "Reaction", "Species"]:
                logger.debug(f"Creating constraint for {label} nodes")
                r = session.run(
                    f"""CREATE CONSTRAINT IF NOT EXISTS
                        ON (n:{label}) ASSERT n.mcId IS UNIQUE;""",
                ).data()
                if r:
                    logger.warning(
                        f"Could not create constraint for {label} nodes: {r}"
                    )
