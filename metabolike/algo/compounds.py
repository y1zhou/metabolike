from difflib import get_close_matches

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
            .assign(value_lowercase=lambda r: r["value"].str.lower())
        )

        return metaid_to_biocyc, res_long

    def search_compound_biocyc_id(self, query: str, **kwargs):
        """Returns the compound biocyc ID if it exists."""
        res = {"query": query, "is_fuzzy_match": False}

        # Direct return if query is a biocyc ID
        if query in self.id_table["biocyc"]:
            res["hits"] = [query]
            return res

        # Try exact matches first
        exact_hits = self.compound_exact_match(query)
        if exact_hits:
            res["hits"] = exact_hits
            return res

        # Try fuzzy search if there's no exact biocyc ID matches
        value_matches = get_close_matches(
            query.lower(), self.cpds["value_lowercase"], **kwargs
        )
        if value_matches:
            hits = [x for v in value_matches for x in self.compound_exact_match(v)]
            res["hits"] = list(set(hits))
            res["is_fuzzy_match"] = True
            return res

        res["hits"] = []
        return res

    def compound_exact_match(self, query: str):
        """Finds biocyc IDs corresponding to the query."""
        res = self.cpds.query("value_lowercase == @query")

        # Return biocyc ID directly if there's an exact match
        if res.shape[0] != 0:
            return list(set(res["biocyc"]))

    def search_compound_compartment(self, query: str):
        """Given a biocyc ID, return compartments it is in."""
        res = self.id_table.loc[self.id_table.biocyc == query, "compartment"]
        return list(set(res.values))

    def search_compound_metaid_in_compartment(
        self, query: str, compartment: str = "cytosol"
    ):
        res = (
            self.id_table.query("biocyc == @query")
            .query("compartment == @compartment")
            .metaId
        )
        return res.tolist()[0]
