import logging
from typing import Any, List, Union

from neo4j import BoltDriver, GraphDatabase, Neo4jDriver, Session, Transaction
from neo4j.work.result import Result

logger = logging.getLogger(__name__)


class BaseDB:
    def __init__(
        self,
        uri: str = "neo4j://localhost:7687",
        neo4j_user: str = "neo4j",
        neo4j_password: str = "neo4j",
    ):
        """setup Neo4j driver.

        Args:
            uri: URI of the Neo4j server. Defaults to ``neo4j://localhost:7687``.
            For more details, see :class:`neo4j.driver.Driver`.
            neo4j_user: Neo4j user. Defaults to "neo4j".
            neo4j_password: Neo4j password. Defaults to "neo4j".

        Attributes:
            driver: :class:`neo4j.Neo4jDriver` or :class:`neo4j.BoltDriver`.
            session: :class:`neo4j.Session`.
        """
        self.driver: Union[Neo4jDriver, BoltDriver] = GraphDatabase.driver(
            uri, auth=(neo4j_user, neo4j_password)
        )
        # logger.debug(self.driver.verify_connectivity())

    def close(self):
        self.session.close()
        self.driver.close()

    def create(self, db_name: str, force: bool = False):
        """Helper function to create a database.

        Args:
            db_name: Name of the database to create.
            force: If True, the database will be dropped if it already exists.
        """
        if not db_name:
            raise ValueError("db_name must be specified")
        if force:
            query = f"""CREATE OR REPLACE DATABASE {db_name}"""
        else:
            query = f"""CREATE DATABASE {db_name} IF NOT EXISTS"""
        with self.driver.session() as ss:
            res = ss.run(query).data()

        if res:
            raise RuntimeError(f"Could not create database {db_name}: {res}")

    def start_session(self, **kwargs):
        self.session: Session = self.driver.session(**kwargs)

    def write(self, cypher: str, **kwargs) -> Result:
        """Helper function to write to the database.

        Args:
            query: Query to write to the database.
            **kwargs: Keyword arguments to pass to :meth:`BaseDB.run`.
        """
        return self.session.write_transaction(lambda tx: tx.run(cypher, **kwargs))

    def read(self, cypher: str, **kwargs) -> Result:
        """Helper function to read from the database.

        Args:
            query: Query to read from the database.
            **kwargs: Keyword arguments to pass to :meth:`BaseDB.run`.
        """
        return self.session.read_transaction(self._read_tx_func, cypher, **kwargs)

    @staticmethod
    def _read_tx_func(tx: Transaction, cypher: str, **kwargs) -> List[Any]:
        """Helper function to run a query.

        Args:
            tx: Neo4j transaction function.
            cypher: Cypher query.
            **kwargs: Keyword arguments to pass to the cypher query.

        Returns:
            The result of Cypher query execution.
            For more details, see :class:`neo4j.work.result.Result`.
        """
        res = tx.run(cypher, **kwargs)
        return [x for x in res]
