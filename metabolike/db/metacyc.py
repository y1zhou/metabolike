import logging
from typing import Any, Dict, List, Tuple

from metabolike.db.base import BaseDB

logger = logging.getLogger(__name__)

NODE_LABELS = [
    "Pathway",
    "Compartment",
    "Reaction",
    "Compound",
    "GeneProduct",
    "RDF",
    "Complex",
    "EntitySet",
    "Citation",
    "Taxa",
]


class MetaDB(BaseDB):
    def use_database(self, db_name: str):
        """Set the name of the database to use."""
        self.db_name = db_name

    def setup_graph_db(self, create_db: bool = True, **kwargs):
        """
        Create Neo4j database and set proper constraints. ``Reaction`` nodes are
        central in the schema.

        Args:
            session: :class:`neo4j.Session` object.
            create_db: If False, does not create the database. This is useful
                for running on neo4j AuraDB when database creation is not allowed.
            **kwargs: Keyword arguments to pass to :meth:`.create`.
        """
        # Create database
        if create_db:
            self.create(self.db_name, **kwargs)

        # Set constraints
        logger.debug("Creating constraint for RDF nodes")
        r = self.write(
            """CREATE CONSTRAINT IF NOT EXISTS
                 ON (r:RDF) ASSERT r.uri IS UNIQUE;""",
        ).data()
        if r:
            logger.warning(f"Could not create constraint for RDF nodes: {r}")

        # Constraints automatically create indexes, so we don't need to
        # create them manually.
        for label in NODE_LABELS:
            logger.debug(f"Creating constraint for node: {label}")
            r = self.write(
                f"""CREATE CONSTRAINT IF NOT EXISTS
                    ON (n:{label}) ASSERT n.mcId IS UNIQUE;""",
            ).data()
            if r:
                logger.warning(f"Could not create constraint for node {label}: {r}")

    def create_node(self, node_label: str, mcid: str, props: Dict[str, str]):
        """Create a node with the given label and properties.

        See :meth:`.Metacyc.sbml_to_graph` for details.

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
        """See :meth:`.Metacyc.sbml_to_graph` for details."""
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
        self, node_label: str, mcid: str, bio_qual: str, props: Dict[str, str]
    ):
        """See :meth:`.Metacyc._add_sbml_rdf_node` for details."""
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
        See :meth:`.Metacyc._add_sbml_gene_product_association_node` for details.

        Link a node to a ``GeneProduct``, ``Complex``, or ``EntitySet`` node.

        Args:
            node_label: Label of the node.
            node_id: ``mcId`` of the node.
            group_id: ``mcId`` of the ``GeneProduct``, ``Complex`` or ``EntitySet``.
            group_label: ``Complex`` or ``EntitySet``.
            edge_type: Type of the edge.
        """
        if group_label not in ["GeneProduct", "Complex", "EntitySet"]:
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

    def get_all_nodes(self, label: str, prop: str) -> List[str]:
        """Fetch an property of nodes with a certain label."""
        if label not in NODE_LABELS:
            raise ValueError(f"Invalid label: {label}")
        # TODO: check for valid properties
        res = self.read(
            f"""
            MATCH (n:{label})
            RETURN DISTINCT n.{prop} AS prop;
            """
        )
        return [n["prop"] for n in res]

    def get_all_compounds(self) -> List[Tuple[str, str]]:
        """Fetch all ``Compound`` nodes.

        All Compound node with RDF have BioCyc IDs.
        The ``META:`` prefix in BioCyc IDs need to be stripped
        """
        res = self.read(
            """
            MATCH (c:Compound)-[:is]->(r:RDF)
            RETURN DISTINCT c.displayName, r.Biocyc;
            """
        )  # TODO: 38 POLYMER nodes don't have BioCyc IDs
        return [(cpd["c.displayName"], cpd["r.Biocyc"][5:]) for cpd in res]

    def link_node_to_citation(
        self, node_label: str, node_display_name: str, citation_id: str
    ):
        """Link a node to a ``Citation`` node.

        Args:
            node_label: Type of the node (Reaction, Pathway, or Compound).
            node_display_name: ``displayName`` of the node.
            citation_id: ``mcId`` of the ``Citation`` node.
        """
        logger.debug(
            f"{node_label} node {node_display_name} has citation {citation_id}"
        )
        self.write(
            f"""
            MATCH (n:{node_label} {{displayName: $dn}})
            MERGE (c:Citation {{mcId: $citation}})
            MERGE (n)-[:hasCitation]->(c)
            """,
            dn=node_display_name,
            citation=citation_id,
        )

    def add_props_to_node(
        self,
        node_label: str,
        node_prop_key: str,
        node_prop_val: str,
        props: Dict[str, Any],
    ):
        """Add properties to a node.

        Args:
            node_label: Label of the node.
            node_prop_key: Property key of the node used to locate the node.
            node_prop_val: Property value of the node used to filter ``node_prop_key``.
            props: Properties to add to the node.
        """
        if node_label not in NODE_LABELS:
            raise ValueError(f"Invalid label: {node_label}")
        # TODO: check for valid properties
        self.write(
            f"""
            MATCH (n:{node_label} {{{node_prop_key}: $val}})
            SET n += $props;
            """,
            val=node_prop_val,
            props=props,
        )

    def add_props_to_compound(
        self, compound_id: str, c_props: Dict[str, Any], rdf_props: Dict[str, Any]
    ):
        self.write(
            """
            MATCH (c:Compound {displayName: $cpd_id})-[:is]->(r:RDF)
            SET c += $c_props, r += $rdf_props;
            """,
            cpd_id=compound_id,
            c_props=c_props,
            rdf_props=rdf_props,
        )

    def link_reaction_to_pathway(self, reaction_id: str, pathway_id: str):
        logger.debug(f"Reaction {reaction_id} is in Pathway {pathway_id}")
        self.write(
            """
            MATCH (r:Reaction {displayName: $reaction})
            MERGE (pw:Pathway {mcId: $pathway})
                ON CREATE SET pw.displayName = $pathway
            MERGE (pw)-[:hasReaction]->(r)
            """,
            reaction=reaction_id,
            pathway=pathway_id,
        )

    def link_pathway_to_pathway(self, pathway_id: str, superpathway_id: str):
        """Link a ``Pathway`` to its super-pathway."""
        logger.debug(f"Pathway {pathway_id} has SuperPathway {superpathway_id}")
        self.write(
            """
            MATCH (pw:Pathway {mcId: $pw_id})
            MERGE (spw:Pathway {mcId: $super_pw})
                ON CREATE SET spw.displayName = $super_pw
            MERGE (spw)-[:hasSubPathway]->(pw)
            """,
            pw_id=pathway_id,
            super_pw=superpathway_id,
        )

    def link_pathway_to_taxa(self, pathway_id: str, tax_id: str, relationship: str):
        """Link a ``Pathway`` to a ``Taxa`` node."""
        if not relationship in {"hasRelatedSpecies", "hasExpectedTaxonRange"}:
            raise ValueError(f"Invalid relationship: {relationship}")
        logger.debug(f"Pathway {pathway_id} {relationship} {tax_id}")
        self.write(
            f"""
            MATCH (pw:Pathway {{mcId: $pw}})
            MERGE (s:Taxa {{mcId: $taxa_id}})
            MERGE (pw)-[:{relationship}]->(s)
            """,
            pw=pathway_id,
            taxa_id=tax_id,
        )

    def link_pathway_to_rate_limiting_step(self, pathway_id: str, reaction_id: str):
        logger.debug(f"Pathway {pathway_id} has rate limiting step {reaction_id}")
        self.write(
            """
            MATCH (pw:Pathway {mcId: $pw}),
                  (r:Reaction {canonical_id: $rxn})
            MERGE (pw)-[l:hasReaction]->(r)
                ON MATCH SET l.isRateLimitingStep = true
            """,
            pw=pathway_id,
            rxn=reaction_id,
        )

    def link_pathway_to_primary_compound(
        self, pathway_id: str, compound_id: str, relationship: str
    ):
        logger.debug(f"Pathway {pathway_id} has primary {relationship} {compound_id}")
        self.write(
            f"""
            MATCH (pw:Pathway {{mcId: $pw}}),
                  (cpd:Compound)-[:is]->(:RDF {{Biocyc: $compound_id}})
            MERGE (pw)-[:hasPrimary{relationship}]->(cpd)
            """,
            pw=pathway_id,
            compound_id=compound_id,
        )

    def link_reaction_to_next_step(self, r1: str, r2: List[str], pathway_id: str):
        """Link a reaction to its next step in a pathway."""
        logger.debug(f"Reaction {r1} has next steps {r2} in {pathway_id}")
        self.write(
            """
            MATCH (r1:Reaction {canonical_id: $r1})
            MATCH (r2:Reaction) WHERE r2.canonical_id IN $r2
            UNWIND r2 AS pred
            MERGE (pred)-[l:isPrecedingEvent]->(r1)
                ON CREATE SET l.hasRelatedPathway = $pw
            """,
            r1=r1,
            r2=r2,
            pw=pathway_id,
        )

    def link_reaction_to_primary_compound(
        self,
        reaction_id: str,
        compound_ids: List[str],
        pathway_id: str,
        side: str,
    ):
        logger.debug(f"Reaction {reaction_id} has primary {side}s {compound_ids}")
        self.write(
            f"""
            MATCH (r:Reaction {{canonical_id: $rxn_id}}),
                  (cpds:Compound)-[:is]->(rdf:RDF)
            WHERE rdf.Biocyc IN $compound_ids
            UNWIND cpds AS cpd
            MATCH (r)-[l:hasLeft|hasRight]->(cpd)
            SET l.isPrimary{side}InPathway = CASE
              WHEN l.isPrimary{side}InPathway IS NULL THEN [$pw]
              WHEN $pw IN l.isPrimary{side}InPathway THEN l.isPrimary{side}InPathway
              ELSE l.isPrimary{side}InPathway + [$pw]
              END;
            """,
            rxn_id=reaction_id,
            compound_ids=compound_ids,
            pw=pathway_id,
        )