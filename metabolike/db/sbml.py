import logging
from typing import Dict, List, Union

from .neo4j import Neo4jClient

logger = logging.getLogger(__name__)


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

    Attributes:
        driver: :class:`neo4j.Neo4jDriver` or :class:`neo4j.BoltDriver`.
        database: str, name of the database to use.
        available_node_labels: tuple of strings indicating the possible node
         labels in the graph.
    """

    def __init__(
        self,
        uri: str = "neo4j://localhost:7687",
        neo4j_user: str = "neo4j",
        neo4j_password: str = "neo4j",
        database: str = "neo4j",
    ):
        super().__init__(uri, neo4j_user, neo4j_password, database)
        self.available_node_labels = (
            "Pathway",
            "Compartment",
            "Reaction",
            "Compound",
            "GeneProduct",
            "RDF",
            "GeneProductComplex",
            "GeneProductSet",
        )

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

        # Set constraints
        logger.debug("Creating constraint for RDF nodes")
        self.write(
            """CREATE CONSTRAINT IF NOT EXISTS ON (r:RDF)
               ASSERT r.uri IS UNIQUE;""",
        )

        # Constraints automatically create indexes, so we don't need to
        # create them manually.
        for label in self.available_node_labels:
            logger.debug(f"Creating constraint for node: {label}")
            self.write(
                f"""CREATE CONSTRAINT IF NOT EXISTS ON (n:{label})
                    ASSERT n.mcId IS UNIQUE;""",
            )

    def create_node(self, node_label: str, mcid: str, props: Dict[str, str]):
        """Create a node with the given label and properties.

        See :meth:`.SBMLParser.sbml_to_graph` for details.

        Args:
            node_label: Label of the node.
            mcid: ``mcId`` of the node.
            props: Properties of the node.
        """
        logger.debug(f"Creating {node_label} node: {mcid}")
        self.write(
            f"""
            MERGE (n:{node_label} {{mcId: $mcId}})
            ON CREATE SET n += $props;
            """,
            mcId=mcid,
            props=props,
        )

    def create_compound_node(self, compartment: str, mcid: str, props: Dict[str, str]):
        """See :meth:`.SBMLParser.sbml_to_graph` for details."""
        logger.debug(f"Creating Compound node: {mcid}")
        self.write(
            """
            MATCH (cpt:Compartment {mcId: $compartment})
            MERGE (c:Compound {mcId: $mcId})-[:hasCompartment]->(cpt)
            ON CREATE SET c += $props;
            """,
            compartment=compartment,
            mcId=mcid,
            props=props,
        )

    def link_node_to_rdf(
        self,
        node_label: str,
        mcid: str,
        bio_qual: str,
        props: Dict[str, Union[str, List[str]]],
    ):
        """See :meth:`.SBMLParser._add_sbml_rdf_node` for details."""
        logger.debug(f"{node_label} node {mcid} {bio_qual} RDF")
        self.write(
            f"""
                MATCH (c:{node_label} {{mcId: $mcId}})
                MERGE (n:RDF)<-[:{bio_qual}]-(c)
                ON CREATE SET n += $props;
                """,
            mcId=mcid,
            props=props,
        )

    def link_reaction_to_compound(
        self,
        reaction_id: str,
        compound_id: str,
        compound_type: str,
        props: Dict[str, str],
    ):
        logger.debug(f"Adding {compound_type} {compound_id} to Reaction {reaction_id}")
        self.write(
            f"""
            MATCH (r:Reaction {{mcId: $reaction}}),
                  (c:Compound {{mcId: $compound}})
            MERGE (r)-[l:has{compound_type}]->(c)
            ON CREATE SET l += $props;
            """,
            reaction=reaction_id,
            compound=compound_id,
            props=props,
        )

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
            node_id: ``mcId`` of the node.
            group_id: ``mcId`` of the ``GeneProduct``, ``GeneProductComplex`` or
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
            MATCH (n:{node_label} {{mcId: $node_id}})
            MERGE (g:{group_label} {{mcId: $group_id}})
            MERGE (g)<-[l:{edge_type}]-(n)
            """,
            node_id=node_id,
            group_id=group_id,
        )
