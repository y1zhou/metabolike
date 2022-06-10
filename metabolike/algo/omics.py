from typing import Dict, Set

from metabolike.db import Neo4jClient


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


def get_table_of_gene_products(db: Neo4jClient, rdf_fields: Dict[str, str]):
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
    for field, col in rdf_fields.items():
        query += f", rdf.{field} AS {col}"
    return db.read(query)
