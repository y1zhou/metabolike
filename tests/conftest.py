from os import environ

import pytest

from metabolike.db.neo4j import Neo4jClient


@pytest.fixture(name="db", scope="session")
def neo4j_connection():
    uri = environ.get("NEO4J_URI")
    username = environ.get("NEO4J_USERNAME", "neo4j")
    password = environ.get("NEO4J_PASSWORD")
    database = environ.get("NEO4J_DATABASE_NAME", "neo4j")

    if not all((uri, username, password)):
        env = "NEO4J_URI, NEO4J_USERNAME, and NEO4J_PASSWORD"
        raise ValueError(f"Environment variables {env} must all be set")

    db = Neo4jClient(
        uri=uri,
        neo4j_user=username,
        neo4j_password=password,
        database=database,
    )
    yield db
    db.close()
