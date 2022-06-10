import logging
from typing import Any, Dict, Iterable, List, Union

from libsbml import FbcModelPlugin, GroupsModelPlugin, Model, Reaction
from metabolike.parser import SBMLParser
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
    MERGE (r)<-[rel:hasRDF]-(gp)
      ON CREATE SET rel.bioQualifier = rdf.bioQual
  );
"""

_gene_product_complex_cypher = """
UNWIND $batch_nodes AS n
  CREATE (gpc:GeneProductComplex:GeneProduct {metaId: n.metaId})
  FOREACH (g IN n.components |
    MERGE (gp:GeneProduct {metaId: g})
    MERGE (gpc)-[:hasComponent]->(gp)
  );
"""

_gene_product_set_cypher = """
UNWIND $batch_nodes AS n
  MERGE (gps:GeneProduct {metaId: n.metaId})
    SET gps:GeneProductSet:GeneProduct
  FOREACH (g IN n.members |
    MERGE (m:GeneProduct {metaId: g})
    MERGE (gps)-[:hasMember]->(m)
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

_reaction_gene_product_cypher = """
UNWIND $batch_nodes AS n
  MERGE (r:Reaction {metaId: n.reaction})
  MERGE (gp:GeneProduct {metaId: n.target})
  MERGE (r)-[:hasGeneProduct]->(gp);
"""

_group_cypher = """
UNWIND $batch_nodes AS n
  MERGE (g:Group {metaId: n.metaId})
    ON CREATE SET g += n.props
  FOREACH (member IN n.members |
    MERGE (m:Reaction {metaId: member})
    MERGE (g)-[:hasGroupMember]->(m)
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

    default_cyphers = {
        "Compartment": _compartment_cypher,
        "Compound": _compound_cypher,
        "GeneProduct": _gene_product_cypher,
        "Reaction": _reaction_cypher,
        "Group": _group_cypher,
        "GeneProductComplex": _gene_product_complex_cypher,
        "GeneProductSet": _gene_product_set_cypher,
        "Reaction-GeneProduct": _reaction_gene_product_cypher,
    }

    def __init__(
        self,
        uri: str = "neo4j://localhost:7687",
        neo4j_user: str = "neo4j",
        neo4j_password: str = "neo4j",
        database: str = "neo4j",
        create_db: bool = True,
        drop_if_exists: bool = False,
        reaction_groups: bool = True,
    ):
        """
        Create Neo4j database.

        Args:
            create_db: If False, does not create the database. This is useful
             for running on neo4j AuraDB when database creation is not allowed.
            drop_if_exists: See :meth:`.create`.
            reaction_groups: Assume all ``group`` nodes are groups of ``Reaction``
                nodes. This assumption greatly speeds up the creation of the
                graph.
        """
        super().__init__(uri, neo4j_user, neo4j_password, database)

        self.available_node_labels = set()
        # Create database
        if create_db:
            self.create(force=drop_if_exists)

        self.reaction_groups = reaction_groups

    def sbml_to_graph(self, parser: SBMLParser):
        """
        Populate Neo4j database with SBML data. The process is as follows:

        #. Parse the SBML file. All parsing errors are logged as warnings.
        #. Create the database and constraints.
        #. Feed the SBML file into the database. This will populate
           ``Compartment``, ``Reaction``, ``Compound``, ``GeneProduct``,
           ``GeneProductSet``, ``GeneProductComplex``, and ``RDF`` nodes.

        Nodes are created for each SBML element using ``MERGE`` statements:
        https://neo4j.com/docs/cypher-manual/current/clauses/merge/#merge-merge-with-on-create
        """
        # Read SBML file
        model: Model = parser.read_sbml(parser.sbml_file).getModel()

        # Compartments
        if cpts := (model.getListOfCompartments()):
            compartments = parser.collect_compartments(cpts)
            self._compartments_to_graph(compartments)

        # Compounds, i.e. metabolites, species
        if cpds := model.getListOfSpecies():
            compounds = parser.collect_compounds(cpds)
            self._compounds_to_graph(compounds)

        # Reactions
        if rxns := (model.getListOfReactions()):
            reactions = parser.collect_reactions(rxns)
            self._reactions_to_graph(reactions)

            # Gene products
            self._gene_products_to_graph(model, parser, rxns)

        # Groups, i.e. related reactions in SBML
        self._groups_to_graph(model, parser)

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
            progress_bar: Show progress bar for slow queries.
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

    def _compartments_to_graph(
        self, compartments: List[Dict[str, Union[str, Dict[str, str]]]]
    ):
        self._set_metaid_constraints("Compartment")
        self.create_nodes(
            "Compartment", compartments, self.default_cyphers["Compartment"]
        )

    def _compounds_to_graph(
        self, compounds: List[Dict[str, Union[str, Dict[str, str]]]]
    ):
        self._set_metaid_constraints("Compound")
        self._check_rdf_nodes(compounds)
        self.create_nodes("Compound", compounds, self.default_cyphers["Compound"])

    def _reactions_to_graph(
        self, reactions: List[Dict[str, Union[str, Dict[str, str]]]]
    ):
        self._set_metaid_constraints("Reaction")
        self._check_rdf_nodes(reactions)
        self.create_nodes(
            "Reaction",
            reactions,
            self.default_cyphers["Reaction"],
            progress_bar=True,
        )

    def _gene_products_to_graph(
        self, model: Model, parser: SBMLParser, reactions: Iterable[Reaction]
    ):
        fbc: FbcModelPlugin = model.getPlugin("fbc")
        if fbc is None:
            logger.warning("No FBC plugin found in SBML file.")
            return

        self._set_metaid_constraints("GeneProduct")
        gene_prods = parser.collect_gene_products(fbc.getListOfGeneProducts())
        self.create_nodes(
            "GeneProduct", gene_prods, self.default_cyphers["GeneProduct"]
        )
        self._check_rdf_nodes(gene_prods)
        self._link_reaction_to_gene_prods(parser, reactions)

    def _link_reaction_to_gene_prods(
        self, parser: SBMLParser, reactions: Iterable[Reaction]
    ):
        (
            reaction_links,
            gene_sets,
            gene_complexes,
        ) = parser.collect_reaction_gene_product_links(reactions)

        # Create complex and set nodes first
        if complex_nodes := (
            [{"metaId": k, "components": list(v)} for k, v in gene_complexes.items()]
        ):
            self._set_metaid_constraints("GeneProductComplex")
            self.create_nodes(
                "GeneProductComplex",
                complex_nodes,
                self.default_cyphers["GeneProductComplex"],
            )

        if set_nodes := (
            [{"metaId": k, "members": list(v)} for k, v in gene_sets.items()]
        ):
            self._set_metaid_constraints("GeneProductSet")
            self.create_nodes(
                "GeneProductSet", set_nodes, self.default_cyphers["GeneProductSet"]
            )

        # Link reactions to these nodes
        reaction_rels = [
            {"reaction": k, "target": v} for k, v in reaction_links.items()
        ]
        self.create_nodes(
            "Reaction-GeneProduct",
            reaction_rels,
            self.default_cyphers["Reaction-GeneProduct"],
        )

    def _groups_to_graph(self, model: Model, parser: SBMLParser):
        groups = model.getPlugin("groups")
        if groups is None:
            logger.warning("No groups plugin found in SBML file.")
        else:
            groups: GroupsModelPlugin
            self._set_metaid_constraints("Group")
            group_nodes = parser.collect_groups(groups.getListOfGroups())

            cypher = self.default_cyphers["Group"]
            if not self.reaction_groups:
                cypher = cypher.replace("m:Reaction", "m")

            self.create_nodes(
                "Group",
                group_nodes,
                cypher,
                batch_size=10,
                progress_bar=True,
            )

    def _set_metaid_constraints(self, label: str):
        """
        Constraints automatically create indexes, so we don't need to
        create them manually.
        """
        if not label:
            raise ValueError("Label must be specified.")

        if label not in self.available_node_labels:
            logger.debug(f"Creating constraint for node: {label}")
            self.write(
                f"""CREATE CONSTRAINT IF NOT EXISTS ON (n:{label})
                    ASSERT n.metaId IS UNIQUE;""",
            )
            self.available_node_labels.add(label)

    def _check_rdf_nodes(self, nodes: Iterable):
        for n in nodes:
            if len(n["rdf"]) > 0:
                self._set_metaid_constraints("RDF")
                break
