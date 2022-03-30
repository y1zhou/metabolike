import logging
from typing import Any, Dict, List, Tuple

import numpy as np
import pandas as pd
from metabolike.utils import generate_gene_reaction_rule

from .sbml import SBMLClient

logger = logging.getLogger(__name__)


# TODO: refactor this with get_high_degree_compound_nodes()
COMMON_COMPOUNDS = ["ATP", "ADP", "H+", "NADH", "NAD+", "H2O", "phosphate"]


class MetacycClient(SBMLClient):
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
            "Citation",
            "Taxa",
        )

    def link_reaction_to_compartment(self, reaction_id: str, compartment_name: str):
        self.write(
            """
            MATCH (r:Reaction {displayName: $reaction}),
                  (c:Compartment {displayName: $compartment})
            MERGE (r)-[:hasCompartment]->(c);
            """,
            reaction=reaction_id,
            compartment=compartment_name,
        )

    def get_all_nodes(self, label: str, prop: str) -> List[str]:
        """Fetch an property of nodes with a certain label."""
        if label not in self.available_node_labels:
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
        """
        Fetch all ``Compound`` nodes. All Compound node with RDF have BioCyc IDs.
        """
        res = self.read(
            """
            MATCH (c:Compound)-[:is]->(r:RDF)
            RETURN DISTINCT c.displayName, r.biocyc;
            """
        )  # TODO: 38 POLYMER nodes don't have BioCyc IDs
        return [(cpd["c.displayName"], cpd["r.biocyc"]) for cpd in res]

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
        if node_label not in self.available_node_labels:
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

    def link_pathway_to_superpathway(self, pathway_id: str, superpathway_id: str):
        """Link a ``Pathway`` to its super-pathway."""
        logger.debug(f"Pathway {pathway_id} has SuperPathway {superpathway_id}")
        self.write(
            """
            MATCH (pw:Pathway {displayName: $pw_id})
            MERGE (spw:Pathway {mcId: $super_pw})
                ON CREATE SET spw.displayName = $super_pw
            MERGE (spw)-[:hasSubPathway]->(pw)
            """,
            pw_id=pathway_id,
            super_pw=superpathway_id,
        )

    def link_pathway_to_pathway(
        self, pw: str, pws: List[str], direction: str, cpd: str
    ):
        """
        Connect a ``Pathway`` to a list of other ``Pathway``s. The ``Compound``
        that is involved in both sides of the connection is specified in the
        relationship.

        We also want to link the ``Reaction`` nodes corresponding to the
        ``Compound`` nodes that are involved in the connection. For example, if
        reaction ``R1`` produces compound ``C``, and reaction ``R2`` consumes
        compound ``C``, we want to link ``R1`` and ``R2`` with a ``isRelatedEvent``
        relationship and add the corresponding pathways as properties.
        """
        if direction == "INCOMING":
            pw_rel_type = "-[l:hasRelatedPathway]->"
            rxn_rel_type = "-[l:isRelatedEvent]->"
            l1, l2 = "hasRight", "hasLeft"
        else:
            pw_rel_type = "<-[l:hasRelatedPathway]-"
            rxn_rel_type = "<-[l:isRelatedEvent]-"
            l1, l2 = "hasLeft", "hasRight"
        self.write(
            f"""
            MATCH (pw1:Pathway {{displayName: $pw_id}})
            MATCH (pw2s:Pathway) WHERE pw2s.displayName IN $pws
            UNWIND pw2s AS pw2
            MERGE (pw1){pw_rel_type}(pw2)
                ON CREATE SET l.hasRelatedCompound = $cpd
            WITH pw1, pw2
            MATCH (pw1)-[:hasReaction]->(r1:Reaction)-[:{l1}]->(:Compound)-[:is]->(:RDF {{biocyc: $cpd}})
            WITH r1, pw2
            MATCH (pw2)-[:hasReaction]->(r2:Reaction)-[:{l2}]->(:Compound)-[:is]->(:RDF {{biocyc: $cpd}})
            MERGE (r1){rxn_rel_type}(r2)
                ON CREATE SET l.fromPathway = $pw_id, l.toPathway = pw2.mcId;
            """,
            pw_id=pw,
            pws=pws,
            cpd=cpd,
        )

    def link_pathway_to_taxa(self, pathway_id: str, tax_id: str, relationship: str):
        """Link a ``Pathway`` to a ``Taxa`` node."""
        if not relationship in {"hasRelatedSpecies", "hasExpectedTaxonRange"}:
            raise ValueError(f"Invalid relationship: {relationship}")
        logger.debug(f"Pathway {pathway_id} {relationship} {tax_id}")
        self.write(
            f"""
            MATCH (pw:Pathway {{displayName: $pw}})
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
            MATCH (pw:Pathway {displayName: $pw}),
                  (r:Reaction {canonicalId: $rxn})
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
            MATCH (cpd:Compound)-[:is]->(:RDF {{Biocyc: $compound_id}}),
                  (pw:Pathway {{displayName: $pw}})-[:hasReaction]->(:Reaction)-[:hasLeft|hasRight]->(cpd)
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
            MATCH (r1:Reaction {canonicalId: $r1})
            MATCH (r2:Reaction) WHERE r2.canonicalId IN $r2
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
            MATCH (r:Reaction {{canonicalId: $rxn_id}}),
                  (cpds:Compound)-[:is]->(rdf:RDF)
            WHERE rdf.biocyc IN $compound_ids
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

    def merge_nodes(
        self, n1_label: str, n2_label: str, n1_attr: str, n2_attr: str, attr_val: str
    ):
        """
        Merge two nodes with the same attribute value. Note that properties is
        hard-coded to be "override", which means all node attributes of ``n2``
        will be used to override the attributes of ``n1``.

        Args:
            n1_label: The label of the first node.
            n2_label: The label of the second node.
            n1_attr: The attribute of the first node for filtering.
            n2_attr: The attribute of the second node for filtering.
            attr_val: The value of the attributes for filtering.
        """
        self.write(
            f"""
            MATCH (n1:{n1_label} {{{n1_attr}: $val}}),
                  (n2:{n2_label} {{{n2_attr}: $val}})
            CALL apoc.refactor.mergeNodes([n1, n2], {{properties: 'override'}})
            YIELD node
            RETURN node;
            """,
            val=attr_val,
        )

    def fix_reaction_direction(self):
        """
        Seems like ``listOfReactants`` in the SBML file are always the true
        reactants in irreversible reactions, no matter if they are under
        ``LEFT`` or ``RIGHT`` in the dat file. Since the reaction nodes were all
        populated using the SBML file, we can say all reactions are either
        reversible or go from left to right.
        """
        self.write(
            """
            MATCH (r:Reaction {reactionDirection: 'right_to_left'})
            SET r.reactionDirection = 'left_to_right'
            """
        )
        self.write(
            """
            MATCH (r:Reaction {reactionDirection: 'physiol_right_to_left'})
            SET r.reactionDirection = 'physiol_left_to_right'
            """
        )

    def delete_bad_reaction_nodes(self):
        """
        Delete reaction nodes with non-canonical mcIds and all their relationships.
        """
        self.write(
            """
            MATCH (r:Reaction) WHERE r.displayName <> r.canonicalId
            DETACH DELETE r;
            """
        )
        self.write(
            """
            MATCH (r:Reaction)
            REMOVE r.canonicalId;
            """
        )

    def set_composite_reaction_labels(self):
        """
        Composite reaction nodes are labeled as ``Pathway`` and ``Reaction``.
        The ``Pathway`` label should be removed and a ``isCompositeReaction``
        property is added to the ``Reaction`` node.
        """
        self.write(
            """
            MATCH (n:Pathway:Reaction)
            REMOVE n:Pathway
            SET n.isCompositeReaction = true;
            """
        )

    def get_all_ec_numbers(self):
        ec = self.read(
            """
            MATCH (r:Reaction)-[:is]->(rdf:RDF)
            WHERE rdf.ecCode IS NOT NULL
            UNWIND rdf.ecCode as ec
            RETURN collect(ec);
            """
        )
        return set(n["ec"] for n in ec)

    def get_view_of_pathway(self, pathway_id: str):
        """
        Get the view of a pathway.

        Args:
            pathway_id: The pathway ID.

        Returns:
            The view of the pathway.
        """
        nodes = self.read(
            """
        MATCH (:Pathway {mcId: $pw_id})-[l:hasReaction]->(r:Reaction)
        WITH r
        MATCH (r)-[:is]->(rrdf:RDF),
              (r)-[:hasGeneProduct|hasComponent|hasMember*]->(gp:GeneProduct)-[:isEncodedBy]->(grdf:RDF)
        WITH r, rrdf.ecCode AS ec,
             COLLECT(gp.displayName) as symbol, COLLECT(grdf.ncbigene) AS ncbi
        RETURN {
          id: id(r), name: r.displayName, ec: ec, gibbs: r.gibbs0,
          direction: r.reactionDirection,
          genes: {symbol: symbol, ncbi: ncbi}
        } AS node;
            """,
            pw_id=pathway_id,
        )

        edges = self.read(
            """
        MATCH (:Pathway {mcId: $pw_id})-[l:hasReaction]->(r:Reaction)
        WITH COLLECT(r) AS nodes
        UNWIND nodes AS r1
        UNWIND nodes AS r2
        WITH * WHERE id(r1) < id(r2)
        MATCH p=(r1)-[:isPrecedingEvent]-(r2)
        UNWIND relationships(p) AS edge
        RETURN {
            id: id(edge), source: id(startNode(edge)), target: id(endNode(edge))
        } AS edge;
            """,
            pw_id=pathway_id,
        )

        return [x["node"] for x in nodes], [x["edge"] for x in edges]

    def get_fba_info_of_pathways(self, pathway_ids: List[str]):
        """
        Retrieve the information of a pathway relevant to flux balance analysis.
        """
        q = self.read(
            """
            MATCH (p:Pathway)-[:hasReaction]->(r:Reaction)-[l:hasLeft|hasRight]->(c:Compound)-[:hasCompartment]->(cpt:Compartment)
            WHERE p.mcId IN $pw_ids
            WITH r, l, c, cpt
            OPTIONAL MATCH (r)-[:is]->(rdf:RDF)
            RETURN
              r.mcId, r.displayName, r.synonyms, r.reactionDirection, r.gibbs0,
              rdf.ecCode,
              type(l), l.stoichiometry,
              c.mcId, c.displayName, c.chemicalFormula, c.gibbs0,
              cpt.mcId;
            """,
            pw_ids=pathway_ids,
        )

        # Use dataframe for easier manipulation
        df = pd.DataFrame(q)
        df.columns = [
            "reaction_id",
            "reaction_name",
            "reaction_synonyms",
            "reaction_direction",
            "reaction_gibbs0",
            "reaction_ec_enzymes",
            "side",
            "stoichiometry",
            "compound_id",
            "compound_name",
            "compound_formula",
            "compound_formation_gibbs0",
            "compound_compartment",
        ]
        all_rxns = df.reaction_id.unique()  # used for getting genes later

        # Fix stoichiometry so that side can be dropped
        df.stoichiometry = np.where(
            df.side == "hasLeft", -df.stoichiometry, df.stoichiometry
        )

        # Fix missing values
        df.reaction_id = df.reaction_name
        df.reaction_name = df.apply(
            lambda x: "|".join(x["reaction_synonyms"])
            if isinstance(x["reaction_synonyms"], list)
            else x["reaction_name"],
            axis=1,
        )
        df.fillna(
            {
                "reaction_gibbs0": pd.NA,
                "compound_formation_gibbs0": pd.NA,
                "reaction_ec_enzymes": pd.NA,
            },
            inplace=True,
        )
        df.drop(columns=["side", "reaction_synonyms"], inplace=True)

        # Get gene associations
        rxn_genes = {}
        for rxn in all_rxns:
            genes = self.get_genes_of_reaction(rxn)
            if genes:
                rxn_genes[rxn] = generate_gene_reaction_rule(genes)

        return df, rxn_genes

    def get_genes_of_reaction(
        self, reaction_id: str
    ) -> List[Tuple[Dict[str, str], Dict[str, str]]]:
        """
        Given a reaction mcId, return the gene products associated with it in
        the form of a list of (source, target) nodes. This is because certain
        reactions are associated with ``GeneProductSet`` or
        ``GeneProductComplex`` nodes, which represent ``OR`` and ``AND``
        relationships, respectively.

        Note that ``GeneProductSet`` could contain nested ``GeneProductComplex``
        nodes.
        """
        res = self.read(
            """
            MATCH (r:Reaction {mcId: $rxn})
            WITH r
            CALL apoc.path.expandConfig(r, {
                labelFilter: "GeneProductSet|GeneProductComplex|/GeneProduct",
                minLevel: 1,
                maxLevel: 3
            })
            YIELD path
            WITH apoc.path.elements(path) AS elements
            UNWIND range(0, size(elements)-2) AS index
            WITH elements, index
            WHERE index % 2 = 0
            RETURN DISTINCT
              {
                label: labels(elements[index])[0],
                mcId: elements[index].mcId,
                name: elements[index].displayName
              } AS source_node,
              {
                label: labels(elements[index+2])[0],
                mcId: elements[index+2].mcId,
                name: elements[index+2].displayName
              } AS target_node;
            """,
            rxn=reaction_id,
        )
        return [(x["source_node"], x["target_node"]) for x in res]

    def get_cpd_view_of_pathway(
        self,
        pathway_id: str,
        ignore_cpds: List[str] = COMMON_COMPOUNDS,
    ):
        """
        Get the view of a pathway with primary compounds linked by reactions.

        Args:
            pathway_id: The pathway mcId or displayName.
            ignore_cpds: A list of common compound displayNames to ignore.
        """
        res = self.read(
            """
            MATCH (:Pathway {displayName: $pw_id})-[:hasReaction]->(r:Reaction),
              (r)-[lr]->(cr:Compound),
              (r)-[lp]->(cp:Compound)
            WHERE (
              ($pw_id IN lr.isPrimaryReactantInPathway) OR
              (r.reactionDirection = 'reversible')
            ) AND (
              ($pw_id IN lp.isPrimaryProductInPathway) OR
              (r.reactionDirection = 'reversible')
            ) AND (
              id(lr) <> id(lp)
            ) AND (
              NOT cr.displayName IN $ignore_cpds
            ) AND (
              NOT cp.displayName IN $ignore_cpds
            )
            RETURN cr, cp,
              apoc.create.vRelationship(cr, r.displayName, properties(r), cp);
            """,
            pw_id=pathway_id,
            ignore_cpds=ignore_cpds,
        )

        return res

    def get_high_degree_compound_nodes(self, degree: int = 60):
        """
        Get all compound nodes with a high degree, i.e. with a lot of edges to
        ``Reaction`` nodes. Including these nodes in many cases would pollute
        the graph with too many outgoing relationships.

        Args:
            degree: The degree threshold.

        Returns:
            A list of mcID of compound nodes.
        """
        return self.read(
            """
            MATCH (c:Compound)<-[:hasLeft|hasRight]-(r:Reaction)
            WITH c.mcId as cpd, COUNT(*) AS cnt
            WHERE cnt >= $degree
            RETURN COLLECT(cpd);
            """,
            degree=degree,
        )

    def get_reaction_route_between_compounds(
        self,
        c1: str,
        c2: str,
        only_pathway_reactions: bool = True,
        ignore_node_mcids: List[str] = [],
        num_routes: int = 2,
        max_hops: int = 10,
    ):
        """
        The function has two modes: one for following pre-defined pathways, and
        one for following any chain of reactions between two compounds.

        Args:
            c1: The first compound mcId.
            c2: The second compound mcId.
            only_pathway_reactions: If True, only follow reactions in pathways.
            ignore_node_mcids: A list of mcIds to ignore.
            num_routes: The number of routes to return.
            max_hops: The maximum number of hops to follow. When argument
             ``only_pathway_reactions`` is False, this is doubled to account
             for the extra hops from reaction nodes to compound nodes.

        Returns:
            A list of possible routes.
        """
        if only_pathway_reactions:
            query = """
            MATCH (n)
            WHERE n.mcId IN $ignore_nodes
            WITH COLLECT(n) AS ns
            MATCH (r:Reaction)-[:hasRight]->(c2:Compound {mcId: $c2})
            WITH c2, ns, COLLECT(r) AS rxns
            MATCH (c1:Compound {mcId: $c1})
            CALL apoc.path.expandConfig(c1, {
                relationshipFilter: "<hasLeft|isPrecedingEvent>|isRelatedEvent>",
                labelFilter: "+Reaction",
                terminatorNodes: rxns,
                blacklistNodes: ns,
                bfs: false,
                limit: $num_routes,
                minLevel: 2,
                maxLevel: $max_hops
            })
            YIELD path
            RETURN c2, path, length(path) AS hops
            ORDER BY hops;
            """
        else:
            query = """
            MATCH (n)
            WHERE n.mcId IN $ignore_nodes
            WITH COLLECT(n) AS ns
            MATCH (c1:Compound {mcId: $c1}),
                  (c2:Compound {mcId: $c2})
            CALL apoc.path.expandConfig(c1, {
                relationshipFilter: "<hasLeft,hasRight>",
                labelFilter: "+Reaction|Compound",
                terminatorNodes: [c2],
                blacklistNodes: ns,
                bfs: true,
                limit: $num_routes,
                minLevel: 2,
                maxLevel: $max_hops
            })
            YIELD path
            RETURN c2, path, length(path) AS hops
            ORDER BY hops;
            """
        res = self.read(
            query,
            ignore_nodes=ignore_node_mcids,
            c1=c1,
            c2=c2,
            num_routes=num_routes,
            max_hops=max_hops * 2,
        )
        return res
