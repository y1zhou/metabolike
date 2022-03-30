import logging
from typing import Any, Callable, Dict, List, Union

from neo4j import BoltDriver, GraphDatabase, Neo4jDriver

logger = logging.getLogger(__name__)


class Neo4jClient:
    """setup Neo4j driver.

    Args:
        uri: URI of the Neo4j server. Defaults to ``neo4j://localhost:7687``.
         For more details, see :class:`neo4j.driver.Driver`.
        neo4j_user: Neo4j user. Defaults to "neo4j".
        neo4j_password: Neo4j password. Defaults to "neo4j".
        database: Name of the database. Defaults to "neo4j".

    Attributes:
        driver: :class:`neo4j.Neo4jDriver` or :class:`neo4j.BoltDriver`.
        database: str, name of the database to use.
    """

    def __init__(
        self,
        uri: str = "neo4j://localhost:7687",
        neo4j_user: str = "neo4j",
        neo4j_password: str = "neo4j",
        database: str = "neo4j",
    ):
        driver = GraphDatabase.driver(uri, auth=(neo4j_user, neo4j_password))
        if not driver:
            raise RuntimeError("Could not connect to Neo4j")
        self.driver: Union[Neo4jDriver, BoltDriver] = driver
        self.database = database
        # logger.debug(self.driver.verify_connectivity())

    def close(self):
        self.driver.close()

    def create(self, force: bool = False):
        """Helper function to create a database.

        Args:
            force: If True, the database will be dropped if it already exists.
        """
        if force:
            query = f"""CREATE OR REPLACE DATABASE {self.database}"""
        else:
            query = f"""CREATE DATABASE {self.database} IF NOT EXISTS"""

        logger.debug(f"Creating database {self.database}")
        with self.driver.session() as ss:
            res = ss.run(query).data()

        if res:
            raise RuntimeError(f"Could not create database {self.database}: {res}")

    def write(self, cypher: str, **kwargs):
        """Helper function to write to the database. Ignores returned output.

        Args:
            query: Query to write to the database.
            **kwargs: Keyword arguments to pass to :meth:`BaseDB.run`.
        """
        with self.driver.session(database=self.database) as ss:
            ss.write_transaction(lambda tx: tx.run(cypher, **kwargs))

    def read(self, cypher: str, **kwargs) -> List[Dict[str, Any]]:
        """Helper function to read from the database. Streams all records in the
        query into a list of dictionaries.

        Args:
            query: Query to read from the database.
            **kwargs: Parameters to pass to the Cypher query.
        """
        with self.driver.session(database=self.database) as ss:
            return ss.read_transaction(lambda tx: tx.run(cypher, **kwargs).data())

    def read_tx(self, tx_func: Callable, **kwargs) -> List[Any]:
        """Helper function to read from the database.

        Args:
            tx_func: A transaction function to run.
            **kwargs: Keyword arguments to pass to ``tx_func`` or parameters for
             the Cypher query.
        """
        with self.driver.session(database=self.database) as ss:
            return ss.read_transaction(tx_func, **kwargs)
