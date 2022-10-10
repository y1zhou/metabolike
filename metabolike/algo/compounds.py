import pandas as pd

from metabolike.db import Neo4jClient


class CompoundMap:
    def __init__(self, db: Neo4jClient):
        self.db = db
        self.id_table, self.cpds = self.get_all_compounds()

    def get_all_compounds(self):
        """
        See https://biocyc.org/PGDBConceptsGuide.shtml#TAG:__tex2page_toc_TAG:__tex2page_sec_4.4
        for detailed descriptions.

        Returns:
            A tuple of two ``pandas.DataFrame``.
              One with columns ``metaId``, ``compartment``, and ``biocyc``;
              One with columns ``biocyc``, ``key_id`` and ``value``.
        """
        cpds = self.db.read(
            """
            MATCH (c:Compound)
            OPTIONAL MATCH (c)-[:hasRDF {bioQualifier: 'is'}]->(rdf:RDF)
            OPTIONAL MATCH (c)-[:hasCompartment]->(cpt:Compartment)
            RETURN c{.metaId, .name, .commonName, .synonyms, .smiles},
                   rdf{.biocyc, .chebi, .keggCompound, .inchikey,
                       .chemspider, .pubchemCompound, .hmdb, .cas,
                       .metabolights, .lipidmaps, .drugbank, .knapsack,
                       .umbbdCompound, .keggGlycan
                   },
                   cpt.metaId AS compartmentId, cpt.commonName AS compartment;
            """
        )  # some don't have the RDF/Compartment fields
        res = [None] * len(cpds)
        for i, cpd in enumerate(cpds):
            entry = cpd["c"]
            if cpd["rdf"] is not None:
                entry = {**entry, **cpd["rdf"]}

            # Prefer compartment commonName, fallback to metaId
            # CCO: Cell Component Ontology; CCO-CYTOSOL is the default location
            if cpd["compartmentId"] is not None:
                entry["compartment"] = cpd["compartmentId"]
            if cpd["compartment"] is not None:
                entry["compartment"] = cpd["compartment"]

            res[i] = entry

        # Convert wide to long table for easier queries
        res = pd.DataFrame(res)
        res = res[~res["biocyc"].isna()]  # POLYMER-xxx, not important
        metaid_to_biocyc = res.filter(["metaId", "biocyc", "compartment"])
        res_long = (
            res.drop(columns=["metaId", "compartment"])
            .melt(id_vars=["biocyc"], var_name="key_id")
            .query("value.notna()")
            .explode("value")
            .drop_duplicates()
        )

        return metaid_to_biocyc, res_long
