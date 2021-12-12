from neo4j import GraphDatabase, Neo4jDriver, Transaction
from neo4j.work.result import Result


def connect(
    host: str = "localhost",
    port: int = 7687,
    neo4j_user: str = "neo4j",
    neo4j_password: str = "neo4j",
) -> Neo4jDriver:
    """setup Neo4j driver.

    Args:
        host: Neo4j host. Defaults to "localhost". Note that there is no need
              to specify the scheme (http:// or bolt://).
        port: Neo4j port. Defaults to 7687.
        neo4j_user: Neo4j user. Defaults to "neo4j".
        neo4j_password: Neo4j password. Defaults to "neo4j".

    Returns:
        Neo4j driver.
    """
    uri = f"{host}:{port}"
    driver = GraphDatabase.neo4j_driver(uri, auth=(neo4j_user, neo4j_password))

    _ = driver.verify_connectivity()
    return driver


def create(driver: Neo4jDriver, db_name: str, force: bool = False):
    """Helper function to create a database.

    Args:
        driver: Neo4j driver created from :func:`connect`.
        db_name: Name of the database to create.
        force: If True, the database will be dropped if it already exists.
    """
    if force:
        query = f"""CREATE OR REPLACE DATABASE {db_name}"""
    else:
        query = f"""CREATE DATABASE {db_name} IF NOT EXISTS"""
    res = driver.session().run(query).data()
    if res:
        raise RuntimeError(f"Could not create database: {db_name} {res}")


def query(tx: Transaction, cypher: str, **kwargs) -> Result:
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
    return res
