from collections.abc import Callable

from neo4j import Neo4jDriver
from pytest_mock import MockerFixture

from metabolike.db.neo4j import Neo4jClient


def test_neo4j_connection(mocker: MockerFixture, db: Neo4jClient):
    assert isinstance(db.driver, Neo4jDriver)
    assert len(db.database) > 0

    # No point checking database creation since it's a demo server
    mocked_session = mocker.patch("neo4j.Session.run")

    db.create()
    mocked_session.assert_called_with(f"CREATE DATABASE {db.database} IF NOT EXISTS")
    db.create(force=True)
    mocked_session.assert_called_with(f"CREATE OR REPLACE DATABASE {db.database}")


def test_read(db: Neo4jClient):
    q = "MATCH (n) RETURN n LIMIT $num_nodes"
    res = db.read(q, num_nodes=10)
    assert len(res) == 10
    assert isinstance(res[0], dict)


def test_read_tx(db: Neo4jClient):
    def tx_func(tx, num_nodes=10):
        resp = tx.run("MATCH (n) RETURN n LIMIT $num_nodes", num_nodes=num_nodes)
        return [x.data() for x in resp]

    res = db.read_tx(tx_func)
    assert len(res) == 10
    assert isinstance(res[0], dict)


def test_write(mocker: MockerFixture, db: Neo4jClient):
    q = "CREATE (n:Movie) SET n.title = $title"
    mocked_session = mocker.patch("neo4j.Session.write_transaction")
    db.write(q, title="The Matrix")

    args, _ = mocked_session.call_args
    assert isinstance(args[0], Callable)
