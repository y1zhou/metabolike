"""
This is NOT usable now. In the future this sub-package will be its own package
that serves the function of proving a web-API for metabolike.
"""

import json
import os

import numpy as np
import pandas as pd
import streamlit as st

from metabolike.algo.compounds import CompoundMap
from metabolike.algo.omics import ReactionGeneMap
from metabolike.algo.routes import find_compound_outflux_routes
from metabolike.db import Neo4jClient

st.title("Route search with metabolike")

# Connect to Neo4j database
neo4j_uri = os.environ.get("NEO4J_URI", "neo4j://localhost:7687")
neo4j_password = os.environ.get("NEO4J_PASSWORD")
neo4j_database = os.environ.get("NEO4J_DATABASE", "neo4j")
con = Neo4jClient(
    uri=neo4j_uri,
    neo4j_password=neo4j_password,
    database=neo4j_database,
)
cpds = CompoundMap(con)
all_cpds = np.concatenate(
    [cpds.id_table["metaId"], cpds.cpds.loc[cpds.cpds["key_id"] == "name", "value"]],
    axis=0,
)
# all_cpds = np.unique(all_cpds)

# Load gene expression data
expression_file = os.environ.get("METABOLIKE_EXP", "./tests/data/transcriptomics.json")
if "gene_exp" not in st.session_state:
    with st.spinner("Loading expression data..."):
        # TODO: use st.selectbox to pick an example dataset
        # TODO: take csv upload with st.file_uploader
        with open(expression_file, "r") as f:
            st.session_state.gene_exp = json.load(f)

rgm = ReactionGeneMap(con, st.session_state.gene_exp)


# Configs for route search
with st.sidebar:
    debug = st.checkbox("Debug mode", value=False)
    st.subheader("Resource limits")
    max_hops = st.slider(
        "Max number of reaction steps", min_value=0, max_value=30, value=15, step=1
    )
    max_num_routes = st.slider(
        "Max number of routes to return", min_value=1, max_value=1000, step=1
    )
    ignore_nodes = st.multiselect(
        "Compounds to ignore",  # TODO: support reaction nodes
        all_cpds,
        default=["H+", "H2O", "chloride", "GDP", "CDP", "ADP", "AMMONIA_c"],
    )

    st.subheader("Route search weights")
    expression_coef = st.slider(
        "Expression level", min_value=0.0, max_value=1.0, value=1.0, step=0.01
    )
    struct_sim_coef = st.slider(
        "End metabolite structure similarity",
        min_value=0.0,
        max_value=1.0,
        value=0.5,
        step=0.01,
    )
    route_len_coef = st.slider(
        "Number of reactions in route",
        min_value=0.0,
        max_value=1.0,
        value=0.0,
        step=0.01,
    )

if debug:
    st.subheader("Expression levels")
    st.dataframe(pd.DataFrame({"logFC": st.session_state.gene_exp}))

# Find compounds
st.subheader("Route search")
cpd_query = st.text_input("Enter compound keyword:")
if cpd_query:
    biocyc_ids = cpds.search_compound_biocyc_id(cpd_query)
    if debug:
        st.write(biocyc_ids)

    cpd_id = st.radio("Compound ID to use:", biocyc_ids["hits"], horizontal=True)
    compartments = cpds.search_compound_compartment(cpd_id)
    compartment = st.radio("Compartment:", compartments, horizontal=True)
    cpd_metaid = cpds.search_compound_metaid_in_compartment(cpd_id, compartment)

    if st.button("Route search"):
        if debug:
            st.write(f"Searching routes starting from {cpd_metaid}")
        routes = find_compound_outflux_routes(
            db=con,
            compound_id=cpd_metaid,
            reaction_gene_expression=rgm,
            drop_nodes=ignore_nodes,
            max_level=max_hops,
            expression_coef=expression_coef,
            structure_similarity_coef=struct_sim_coef,
            route_length_coef=route_len_coef,
        )

        if routes:
            st.write(
                f"Found {len(routes)} routes in total, returning top {max_num_routes}"
            )
            for r in routes:
                r["genes"] = "|".join(list(r["genes"]))

            st.dataframe(
                pd.DataFrame(routes)
                .assign(
                    chain=lambda df: df["route"].apply(
                        lambda r: " -> ".join(
                            x["name"] for x in r if x["nodeType"] == "Compound"
                        )
                    )
                )
                .drop(columns="route")
                .drop_duplicates()
                .head(max_num_routes),
                use_container_width=True,
            )

con.close()
