from typing import Dict, List, Set, Tuple

import numpy as np
import pandas as pd
from metabolike.db import Neo4jClient
from metabolike.utils import generate_gene_reaction_rule

# TODO: refactor this with get_high_degree_compound_nodes()
COMMON_COMPOUNDS = ["ATP", "ADP", "H+", "NADH", "NAD+", "H2O", "phosphate"]


def get_all_ec_numbers(db: Neo4jClient) -> Set[str]:
    ec = db.read(
        """
        MATCH (r:Reaction)-[:hasRDF {bioQualifier: 'is'}]->(rdf:RDF)
        WHERE rdf.ecCode IS NOT NULL
        UNWIND rdf.ecCode as ec
        RETURN collect(ec);
        """
    )
    return {n["ec"] for n in ec}


def get_view_of_pathway(db: Neo4jClient, pathway_id: str):
    """
    Get the view of a pathway.

    Args:
        pathway_id: The pathway ID.

    Returns:
        The view of the pathway.
    """
    nodes = db.read(
        """
    MATCH (:Pathway {metaId: $pw_id})-[l:hasReaction]->(r:Reaction)
    WITH r
    MATCH (r)-[:hasRDF {bioQualifier: 'is'}]->(rrdf:RDF),
            (r)-[:hasGeneProduct|hasComponent|hasMember*]->(gp:GeneProduct)-[:hasRDF {bioQualifier: 'isEncodedBy'}]->(grdf:RDF)
    WITH r, rrdf.ecCode AS ec,
            COLLECT(gp.name) as symbol, COLLECT(grdf.ncbigene) AS ncbi
    RETURN {
        id: id(r), name: r.name, ec: ec, gibbs: r.gibbs0,
        direction: r.reactionDirection,
        genes: {symbol: symbol, ncbi: ncbi}
    } AS node;
        """,
        pw_id=pathway_id,
    )

    edges = db.read(
        """
    MATCH (:Pathway {metaId: $pw_id})-[l:hasReaction]->(r:Reaction)
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


def get_fba_info_of_pathways(db: Neo4jClient, pathway_ids: List[str]):
    """
    Retrieve the information of a pathway relevant to flux balance analysis.
    """
    q = db.read(
        """
    MATCH (p:Pathway)-[:hasReaction]->(r:Reaction)-[l:hasLeft|hasRight]->(c:Compound)-[:hasCompartment]->(cpt:Compartment)
    WHERE p.metaId IN $pw_ids
    WITH r, l, c, cpt
    OPTIONAL MATCH (r)-[:hasRDF {bioQualifier: 'is'}]->(rdf:RDF)
    RETURN
        r.metaId, r.name, r.synonyms, r.reactionDirection, r.gibbs0,
        rdf.ecCode,
        type(l), l.stoichiometry,
        c.metaId, c.name, c.chemicalFormula, c.gibbs0,
        cpt.metaId;
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
        if genes := get_genes_of_reaction(db, rxn):
            rxn_genes[rxn] = generate_gene_reaction_rule(genes)

    return df, rxn_genes


def get_genes_of_reaction(
    db: Neo4jClient, reaction_id: str
) -> List[Tuple[Dict[str, str], Dict[str, str]]]:
    """
    Given a reaction metaId, return the gene products associated with it in
    the form of a list of (source, target) nodes. This is because certain
    reactions are associated with ``GeneProductSet`` or
    ``GeneProductComplex`` nodes, which represent ``OR`` and ``AND``
    relationships, respectively.

    Note that ``GeneProductSet`` could contain nested ``GeneProductComplex``
    nodes.
    """
    res = db.read(
        """
        MATCH (r:Reaction {metaId: $rxn})
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
            metaId: elements[index].metaId,
            name: elements[index].name
            } AS source_node,
            {
            label: labels(elements[index+2])[0],
            metaId: elements[index+2].metaId,
            name: elements[index+2].name
            } AS target_node;
        """,
        rxn=reaction_id,
    )
    return [(x["source_node"], x["target_node"]) for x in res]


def get_cpd_view_of_pathway(
    db: Neo4jClient,
    pathway_id: str,
    ignore_cpds: List[str] = COMMON_COMPOUNDS,
):
    """
    Get the view of a pathway with primary compounds linked by reactions.

    Args:
        pathway_id: The pathway metaId or name.
        ignore_cpds: A list of common compound names to ignore.
    """
    return db.read(
        """
        MATCH (:Pathway {name: $pw_id})-[:hasReaction]->(r:Reaction),
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
            NOT cr.name IN $ignore_cpds
        ) AND (
            NOT cp.name IN $ignore_cpds
        )
        RETURN cr, cp,
            apoc.create.vRelationship(cr, r.name, properties(r), cp);
        """,
        pw_id=pathway_id,
        ignore_cpds=ignore_cpds,
    )


def get_high_degree_compound_nodes(db: Neo4jClient, degree: int = 60):
    """
    Get all compound nodes with a high degree, i.e. with a lot of edges to
    ``Reaction`` nodes. Including these nodes in many cases would pollute
    the graph with too many outgoing relationships.

    Args:
        degree: The degree threshold.

    Returns:
        A list of mcID of compound nodes.
    """
    return db.read(
        """
        MATCH (c:Compound)<-[:hasLeft|hasRight]-(r:Reaction)
        WITH c.metaId as cpd, COUNT(*) AS cnt
        WHERE cnt >= $degree
        RETURN COLLECT(cpd);
        """,
        degree=degree,
    )


def get_reaction_route_between_compounds(
    db: Neo4jClient,
    c1: str,
    c2: str,
    only_pathway_reactions: bool = True,
    ignore_node_metaids: List[str] = [],
    num_routes: int = 2,
    max_hops: int = 10,
):
    """
    The function has two modes: one for following pre-defined pathways, and
    one for following any chain of reactions between two compounds.

    Args:
        c1: The first compound metaId.
        c2: The second compound metaId.
        only_pathway_reactions: If True, only follow reactions in pathways.
        ignore_node_metaids: A list of metaIds to ignore.
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
        WHERE n.metaId IN $ignore_nodes
        WITH COLLECT(n) AS ns
        MATCH (r:Reaction)-[:hasRight]->(c2:Compound {metaId: $c2})
        WITH c2, ns, COLLECT(r) AS rxns
        MATCH (c1:Compound {metaId: $c1})
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
        WHERE n.metaId IN $ignore_nodes
        WITH COLLECT(n) AS ns
        MATCH (c1:Compound {metaId: $c1}),
                (c2:Compound {metaId: $c2})
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
    return db.read(
        query,
        ignore_nodes=ignore_node_metaids,
        c1=c1,
        c2=c2,
        num_routes=num_routes,
        max_hops=max_hops * 2,
    )


def find_compound_sink(
    db: Neo4jClient,
    compound_id: str,
    expressed_genes: List[str],
    max_hops: int = 10,
):
    """
    Find how a compound is consumed. Genes are either expressed or unexpressed,
    so this is a rather binary view of the metabolic network.

    Args:
        db: Database client.
        compound_id: The metaId of the compound of interest.
        expressed_genes: A list of genes that are considered expressed.
        max_hops: The maximum number of downstream compounds before stopping.
    """
    return db.read(
        """
        // List of considered genes
        MATCH (r:Reaction)-[:hasGeneProduct|hasMember|hasComponent*]->(gp:geneProduct)
        WHERE gp.id IN $genes
        WITH COLLECT(r) as valid_rxns
        MATCH (c:Compound {metaId: $cpd_id})<-[l:hasLeft|hasRight]-(r:Reaction)
        // keep reversible reactions and reactions that consume the cpd
        WHERE NOT (
            type(l) = 'hasRight' AND
            r.reactionDirection CONTAINS 'left_to_right'
        )
        WITH r, c
        RETURN r, c
        """,
        cpd_id=compound_id,
        genes=expressed_genes,
        max_hops=max_hops,
    )


def get_all_reactions_related_to_compound(
    db: Neo4jClient, cpd_name: str
) -> pd.DataFrame:
    q = db.read(
        """
        MATCH (r:Reaction)-[l:hasLeft|hasRight]->(c:Compound {name: $cpd_name})
        WITH r, type(l) AS side, l.stoichiometry AS stoichiometry
        MATCH (r)-[:hasGeneProduct|hasComponent|hasMember*]->(g:GeneProduct)-[:hasRDF {bioQualifier: 'isEncodedBy'}]->(rdf:RDF)
        WITH r, side, stoichiometry,
          COLLECT(g.name) AS gene_symbols, COLLECT(COALESCE(rdf.ncbigene, 'MISSING')) AS entrez_ids
        RETURN
          r.name AS reaction_id,
          r.reactionDirection AS reaction_direction,
          r.types AS reaction_types,
          side AS compound_on_side,
          stoichiometry,
          gene_symbols,
          entrez_ids;
        """,
        cpd_name=cpd_name,
    )

    # Clean up missing values and annotate transport reactions
    res = pd.DataFrame(list(q))
    res["reaction_direction"] = res["reaction_direction"].fillna("unknown")
    res["is_transport_reaction"] = res["reaction_types"].map(
        lambda v: any(x.startswith("TR-") for x in v)
    )
    res.replace(
        {"compound_on_side": {"hasLeft": "left", "hasRight": "right"}}, inplace=True
    )

    res.drop(columns=["reaction_types"], inplace=True)
    res.sort_values(
        ["is_transport_reaction", "reaction_direction", "compound_on_side"],
        ignore_index=True,
        inplace=True,
    )

    return res


def export_pathways(db: Neo4jClient) -> pd.DataFrame:
    q = db.read(
        """
        MATCH (p:Pathway)-[:hasReaction]->(:Reaction)-[:hasGeneProduct|hasComponent|hasMember*]->(g:GeneProduct)
        RETURN p.metaId AS pathway_id,
          p.commonName AS pathway_name,
          p.types AS pathway_type,
          p.synonyms AS pathway_synonyms,
          COLLECT(g.name) AS gene_symbols;
        """
    )
    return pd.DataFrame(list(q))
