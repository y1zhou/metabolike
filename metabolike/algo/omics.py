import logging
from typing import Dict, Iterable, Sequence, Set

import pandas as pd

from metabolike.db import Neo4jClient

logger = logging.getLogger(__name__)


def get_all_gene_products(db: Neo4jClient) -> Set[str]:
    ec = db.read(
        """
        MATCH (r:Reaction)-[:hasRDF {bioQualifier: 'is'}]->(rdf:RDF)
        WHERE rdf.ecCode IS NOT NULL
        UNWIND rdf.ecCode as ec
        RETURN collect(ec);
        """
    )
    return {n["ec"] for n in ec}


def get_table_of_gene_products(db: Neo4jClient, rdf_fields: Dict[str, str] = None):
    """
    Retrieves all reaction-associated gene products.
    Args:
        db: A Neo4j client connected to the graph database.
        rdf_fields: Properties of the RDF nodes, and the desired output name.
          For example, {"ncbigene": "entrez"} would extract the ``ncbigene``
          property from the RDF nodes and output it in the ``entrez`` column.

    Returns:
        A list of entries with the reaction ID, the metaId of the gene, the
          gene symbol, and other fields from RDF nodes whenever available.
    """
    query = """
        MATCH (r:Reaction)-[:hasGeneProduct|hasMember|hasComponent*]->(gp:GeneProduct)
        WHERE NOT gp:GeneProductSet OR gp:GeneProductComplex
        WITH r.metaId AS rxn_id, gp
        OPTIONAL MATCH (gp)-[:hasRDF {bioQualifier: 'isEncodedBy'}]->(rdf:RDF)
        RETURN rxn_id, gp.metaId AS gene_id, gp.name AS gene
    """
    if rdf_fields:
        for field, col in rdf_fields.items():
            query += f", rdf.{field} AS {col}"
    return db.read(query)


class ReactionGeneMap:
    def __init__(
        self,
        database_connection: Neo4jClient,
        gene_ids: Sequence[str],
        expression_levels: Sequence[float],
        gene_groups: pd.DataFrame | None = None,
        gene_id_map: pd.DataFrame | None = None,
        gene_set_reduce_func: callable = max,
        gene_complex_reduce_func: callable = min,
    ):
        """

        Attributes:
            database_connection: Connection to the Neo4j database.
            rxn_exp: expression levels of each reaction.
            gene_ids: A list-like of strings contains gene identifiers.
            expression_levels: gene expression data.
            gene_groups: parent-child mapping of gene product sets and complexes.
            gene_id_map: mapping between gene metaId and gene name.
        """

        self.db = database_connection
        self.gene_exp = {
            gene_id: exp_level
            for gene_id, exp_level in zip(gene_ids, expression_levels)
        }
        self.rxn_exp: dict[str, float] = {}

        if gene_groups is None:
            self.gene_groups = self._get_gene_groups()
        else:
            self.gene_groups = gene_groups

        if gene_id_map is None:
            self.gene_id_map = self._get_gene_id_map()
        else:
            self.gene_id_map = gene_id_map

        self._valid_gene_groups = set(self.gene_groups["group_id"].unique())
        self._set_func = gene_set_reduce_func
        self._complex_func = gene_complex_reduce_func

        # Setup reaction expression levels
        rxn_genes = self._get_top_level_rxn_gene_mapping()
        self.reaction_gp_mapping = {}

        for rxn, gene in rxn_genes.items():
            # GeneProductSet or GeneProductComplex
            if gene in self._valid_gene_groups:
                self.reaction_gp_mapping[rxn] = self.get_all_genes_in_group(gene)
                self.rxn_exp[rxn] = self.get_gene_expression(gene)
            elif gene_name := self.gene_id_map.get(gene):  # not None or ""
                self.reaction_gp_mapping[rxn] = {gene_name}
                self.rxn_exp[rxn] = self.get_gene_expression(gene_name)

    def get_all_genes_in_group(self, gene: str):
        if gene in self.gene_exp:
            return {gene}
        elif gene in self._valid_gene_groups:
            res = set()
            genes = self.gene_groups.loc[self.gene_groups["group_id"] == gene, :]
            gene_ids = genes["members"].values
            for g in gene_ids:
                if gene_name := self.gene_id_map.get(g):
                    if maybe_gene_name := self.get_all_genes_in_group(gene_name):
                        res.union(maybe_gene_name)

            return res

    def get_gene_expression(self, gene_name: str):
        if (gene_exp := self.gene_exp.get(gene_name)) is not None:
            return gene_exp
        elif gene_name in self._valid_gene_groups:
            return self.calc_reaction_gene_group_expression(gene_name)
        # else:
        #     # raise ValueError(f"Unknown gene product {gene_name}")
        #     logger.warning(f"Unknown gene product {gene_name}")

    def calc_reaction_gene_group_expression(self, gene_group_id: str):
        """
        Calculate the expression level of a given gene product set or
        complex node.
        """
        # TODO: make use of self.get_all_genes_in_group()
        genes = self.gene_groups.loc[self.gene_groups["group_id"] == gene_group_id, :]
        gene_ids = genes["members"].values
        group_type = genes["group_type"].values[0]

        res = [
            self.get_gene_expression(gene_name)
            for gene_id in gene_ids
            if (gene_name := self.gene_id_map.get(gene_id))
        ]
        res = [x for x in res if x is not None]
        if not res:
            return 0.0

        if group_type == "Set":
            return self._set_func(res)
        elif group_type == "Complex":
            return self._complex_func(res)
        else:
            raise ValueError(f"Unknown gene group type: {group_type}")

    def get_route_expression(self, rxn_ids: Iterable[str]) -> float:
        """
        Sums up total expression levels of the reaction route, and divide
        the value by the number of reactions.

        Args:
            rxn_ids: The metaId fields of the reactions.

        Returns:
            The averaged expression level.
        """
        res = 0.0
        for rxn_id in rxn_ids:
            if (gene_exp := self.rxn_exp.get(rxn_id)) is not None:
                res += gene_exp

        return res / len(rxn_ids)

    def _get_top_level_rxn_gene_mapping(self):
        """
        Get unique reaction -> gene product mappings. The "gene" column could
        be GeneProduct, GeneProductSet, or GeneProductComplex nodes.
        """
        rxn2gene = self.db.read(
            """
            MATCH (r:Reaction)-[:hasGeneProduct]->(g)
            WHERE g:GeneProduct OR g:GeneProductSet OR g:GeneProductComplex
            RETURN r.metaId AS rxn, g.metaId AS gene_id;
            """
        )
        res = {n["rxn"]: n["gene_id"] for n in rxn2gene}
        return res

    def _get_gene_groups(self):
        group2child = self.db.read(
            """
            MATCH (gps:GeneProductSet)-[:hasMember]->(gp)
            RETURN gps.metaId AS group_id, "Set" AS group_type, COLLECT(gp.metaId) AS members
            UNION ALL
            MATCH (gpc:GeneProductComplex)-[:hasComponent]->(gp)
            RETURN gpc.metaId AS group_id, "Complex" AS group_type, COLLECT(gp.metaId) AS members
            """
        )
        res = pd.DataFrame([x for x in group2child])
        return res.explode("members").reset_index(drop=True)

    def _get_gene_id_map(self):
        genes = self.db.read(
            """
            MATCH (gp:GeneProduct)
            WHERE NOT (gp:GeneProductSet OR gp:GeneProductComplex)
            RETURN gp.metaId AS gene_id, gp.name AS gene_name;
            """
        )
        return {n["gene_id"]: n["gene_name"] for n in genes}
