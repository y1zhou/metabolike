import logging
from collections.abc import Callable
from typing import Any, Union

from neo4j import BoltDriver, GraphDatabase, Neo4jDriver

logger = logging.getLogger(__name__)


class Neo4jClient:
    """setup Neo4j driver.

    Args:
        uri: URI of the Neo4j server. Defaults to `neo4j://localhost:7687`.
            For more details, see [neo4j.driver.Driver](https://neo4j.com/docs/api/python-driver/current/api.html).
        neo4j_user: Neo4j username.
        neo4j_password: Neo4j password.
        database: Name of the database.

    Attributes:
        driver: `neo4j.Neo4jDriver` or `neo4j.BoltDriver`.
        database: name of the database to use.
    """

    def __init__(  # nosec B107 - default password is okay here
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
        self.database: str = database

        _connected = self.driver.get_server_info()
        logger.debug(f"Connected to {_connected.address}")

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
            ss.run(query)

    def write(self, cypher: str, **kwargs):
        """Helper function to write to the database. Ignores returned output.

        Args:
            cypher: Query to write to the database.
            **kwargs: Keyword arguments to pass to :meth:`BaseDB.run`.
        """
        with self.driver.session(database=self.database) as ss:
            ss.execute_write(lambda tx: tx.run(cypher, **kwargs))

    def read(self, cypher: str, **kwargs) -> list[dict[str, Any]]:
        """Helper function to read from the database. Streams all records in the query into a list
        of dictionaries.

        Args:
            cypher: Query to read from the database.
            **kwargs: Parameters to pass to the Cypher query.
        """
        with self.driver.session(database=self.database) as ss:
            return ss.execute_read(lambda tx: tx.run(cypher, **kwargs).data())

    def read_tx(self, tx_func: Callable, **kwargs) -> list[Any]:
        """Helper function to read from the database.

        Args:
            tx_func: A transaction function to run.
            **kwargs: Keyword arguments to pass to ``tx_func`` or parameters for
                the Cypher query.
        """
        with self.driver.session(database=self.database) as ss:
            return ss.execute_read(tx_func, **kwargs)
