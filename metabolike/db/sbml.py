import logging
from typing import Any, Dict, List

from metabolike.utils import chunk
from tqdm import tqdm

from .neo4j import Neo4jClient

logger = logging.getLogger(__name__)

_compartment_cypher = """
UNWIND $batch_nodes AS n
  MERGE (x:Compartment {metaId: n.metaId})
    ON CREATE SET x += n.props;
"""

_compound_cypher = """
UNWIND $batch_nodes AS n
  MATCH (cpt:Compartment {metaId: n.compartment})
  MERGE (c:Compound {metaId: n.metaId})-[:hasCompartment]->(cpt)
    ON CREATE SET c += n.props
  FOREACH (rdf IN n.rdf |
    CREATE (r:RDF)
    SET r = rdf.rdf
    MERGE (r)<-[rel:hasRDF]-(c)
      ON CREATE SET rel.bioQualifier = rdf.bioQual
  );
"""

_gene_product_cypher = """
UNWIND $batch_nodes AS n
  MERGE (gp:GeneProduct {metaId: n.metaId})
    ON CREATE SET gp += n.props
  FOREACH (rdf IN n.rdf |
    CREATE (r:RDF)
    SET r = rdf.rdf
    MERGE (r)<-[rel:hasRDF]-(c)
      ON CREATE SET rel.bioQualifier = rdf.bioQual
  );
"""

_reaction_cypher = """
UNWIND $batch_nodes AS n
  MERGE (r:Reaction {metaId: n.metaId})
    ON CREATE SET r += n.props
  FOREACH (rdf IN n.rdf |
    CREATE (x:RDF)
    SET x = rdf.rdf
    MERGE (x)<-[rel:hasRDF]-(r)
      ON CREATE SET rel.bioQualifier = rdf.bioQual
  )
  FOREACH (reactant IN n.reactants |
    MERGE (c:Compound {metaId: reactant.cpdId})  // can't MATCH here
    MERGE (r)-[rel:hasLeft]->(c)
      ON CREATE SET rel = reactant.props
  )
  FOREACH (product IN n.products |
    MERGE (c:Compound {metaId: product.cpdId})
    MERGE (r)-[rel:hasRight]->(c)
      ON CREATE SET rel = product.props
  );
"""


class SBMLClient(Neo4jClient):
    """In addition to the Neo4j driver, this class also includes a set of
    helper methods for creating nodes and relationships in the graph with data
    from the SBML file.

    Args:
        uri: URI of the Neo4j server. Defaults to ``neo4j://localhost:7687``.
         For more details, see :class:`neo4j.driver.Driver`.
        neo4j_user: Neo4j user. Defaults to ``neo4j``.
        neo4j_password: Neo4j password. Defaults to ``neo4j``.
        database: Name of the database. Defaults to ``neo4j``.
        create_db: Whether to create the database. See :meth:`.setup_graph_db`.
        drop_if_exists: Whether to drop the database if it already exists.

    Attributes:
        driver: :class:`neo4j.Neo4jDriver` or :class:`neo4j.BoltDriver`.
        database: str, name of the database to use.
        available_node_labels: tuple of strings indicating the possible node
         labels in the graph.
    """

    available_node_labels = (
        "Pathway",
        "Compartment",
        "Reaction",
        "Compound",
        "GeneProduct",
        "RDF",
        "GeneProductComplex",
        "GeneProductSet",
    )

    default_cyphers = {
        "Compartment": _compartment_cypher,
        "Compound": _compound_cypher,
        "GeneProduct": _gene_product_cypher,
        "Reaction": _reaction_cypher,
    }

    def __init__(
        self,
        uri: str = "neo4j://localhost:7687",
        neo4j_user: str = "neo4j",
        neo4j_password: str = "neo4j",
        database: str = "neo4j",
        create_db: bool = True,
        drop_if_exists: bool = False,
    ):
        super().__init__(uri, neo4j_user, neo4j_password, database)
        self.setup_graph_db(create_db=create_db, drop_if_exists=drop_if_exists)

    def setup_graph_db(self, create_db: bool = True, drop_if_exists: bool = False):
        """
        Create Neo4j database and set proper constraints. ``Reaction`` nodes are
        central in the schema.

        Args:
            session: :class:`neo4j.Session` object.
            create_db: If False, does not create the database. This is useful
             for running on neo4j AuraDB when database creation is not allowed.
            drop_if_exists: See :meth:`.create`.
        """
        # Create database
        if create_db:
            self.create(force=drop_if_exists)

        # Constraints automatically create indexes, so we don't need to
        # create them manually.
        for label in self.available_node_labels:
            logger.debug(f"Creating constraint for node: {label}")
            self.write(
                f"""CREATE CONSTRAINT IF NOT EXISTS ON (n:{label})
                    ASSERT n.metaId IS UNIQUE;""",
            )

    def create_nodes(
        self,
        desc: str,
        nodes: List[Dict[str, Any]],
        query: str,
        batch_size: int = 1000,
        progress_bar: bool = False,
    ):
        """Create nodes in batches with the given label and properties.

        For ``Compartment`` nodes, simply create them with given properties.

        Each ``compound`` node is linked to its ``Compartment`` node. If it
        has related ``RDF`` nodes, these are also created and linked to the
        ``Compound`` node.

        ``GeneProduct`` nodes don't have relationships to ``Compartment``
        nodes, but they are linked to corresponding ``RDF`` nodes.

        Args:
            desc: Label of the node in log and progress bar.
            nodes: List of properties of the nodes.
            query: Cypher query to create the nodes.
            batch_size: Number of nodes to create in each batch.
        """
        logger.info(f"Creating {desc} nodes")

        if progress_bar:
            it = tqdm(
                chunk(nodes, batch_size),
                desc=desc,
                total=len(nodes) // batch_size,
            )
        else:
            it = chunk(nodes, batch_size)

        for batch in it:
            self.write(query, batch_nodes=batch)

        logger.info(f"Created {len(nodes)} {desc} nodes")

    def link_node_to_gene(
        self,
        node_label: str,
        node_id: str,
        group_id: str,
        group_label: str,
        edge_type: str,
    ):
        """
        See :meth:`.SBMLParser._add_sbml_gene_product_association_node` for
        details.

        Link a node to a ``GeneProduct``, ``GeneProductComplex``, or
        ``GeneProductSet`` node.

        Args:
            node_label: Label of the node.
            node_id: ``metaId`` of the node.
            group_id: ``metaId`` of the ``GeneProduct``, ``GeneProductComplex`` or
             ``GeneProductSet``.
            group_label: ``GeneProductComplex`` or ``GeneProductSet``.
            edge_type: Type of the edge.
        """
        if group_label not in {"GeneProduct", "GeneProductComplex", "GeneProductSet"}:
            raise ValueError(f"Invalid group_type: {group_label}")
        logger.debug(
            f"{node_label} node {node_id} {edge_type} {group_label} {group_id}"
        )
        self.write(
            f"""
            MATCH (n:{node_label} {{metaId: $node_id}})
            MERGE (g:{group_label} {{metaId: $group_id}})
            MERGE (g)<-[l:{edge_type}]-(n)
            """,
            node_id=node_id,
            group_id=group_id,
        )
