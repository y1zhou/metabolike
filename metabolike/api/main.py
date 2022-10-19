"""
This is NOT usable now. In the future this sub-package will be its own package
that serves the function of proving a streamlit app for route searching.
"""

import numpy as np
import pandas as pd
import streamlit as st

from metabolike.algo.omics import ReactionGeneMap
from metabolike.algo.routes import find_compound_outflux_routes
from metabolike.api.data import (
    TCGA,
    download_df,
    download_tcga_data,
    get_compound_map,
    neo4j_db,
)

st.title("Route search with metabolike")

# Connect to Neo4j database
con = neo4j_db()

# Load dict of all compounds in graph
cpds = get_compound_map(con)
all_cpds = np.concatenate(
    [cpds.cpds.loc[cpds.cpds["key_id"] == "name", "value"], cpds.id_table["metaId"]],
    axis=0,
)

# Configs for route search
with st.sidebar:
    st.subheader("Resource limits")
    max_hops = st.slider(
        "Max number of reaction steps", min_value=0, max_value=30, value=15, step=1
    )
    max_num_routes = st.slider(
        "Max number of routes to return", min_value=1, max_value=200, value=20, step=1
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
        value=0.0,
        step=0.01,
    )
    route_len_coef = st.slider(
        "Number of reactions in route",
        min_value=0.0,
        max_value=1.0,
        value=0.0,
        step=0.01,
    )

    st.subheader("Route search options")
    route_mode = st.radio(
        "Top routes", options=("Highest scores", "Lowest scores"), horizontal=True
    )
    debug = st.checkbox("Debug mode", value=False)


# Load expression dataset
st.subheader("Load expression data")
col1, col2 = st.columns([2, 3], gap="large")

if "exp_loaded" not in st.session_state:
    st.session_state.exp_loaded = False
with col1:
    tcga = st.selectbox(
        "You can test the route search function with a TCGA example dataset:",
        options=TCGA,
    )

with col2:
    uploaded_file = st.file_uploader(
        "Alternatively, upload your own dataset:", type="csv"
    )


if st.button("Load dataset"):
    if uploaded_file is not None:
        with st.spinner("Loading uploaded dataset..."):
            st.session_state.gene_exp = pd.read_csv(uploaded_file)
            st.session_state.gene_id_col = st.selectbox(
                "Column for gene IDs:",
                options=st.session_state.gene_exp.columns,
                index=0,
            )
            st.session_state.score_col = st.selectbox(
                "Column for scores:", options=st.session_state.gene_exp.columns, index=1
            )
    else:
        with st.spinner("Loading TCGA dataset..."):
            st.session_state.gene_exp = download_tcga_data(tcga)
            st.session_state.gene_id_col = "Symbol"
            st.session_state.score_col = "Score"
    st.session_state.exp_loaded = True

if st.session_state.exp_loaded:
    if debug:
        st.write("Gene ID column: ", st.session_state.gene_id_col)
        st.write("Gene scores column: ", st.session_state.score_col)
    gene_ids = st.session_state.gene_exp[st.session_state.gene_id_col].to_list()
    exp_levels = st.session_state.gene_exp[st.session_state.score_col].to_list()
    rgm = ReactionGeneMap(con, gene_ids, exp_levels)
    st.success("Expression data loaded!")

    if debug:
        st.dataframe(st.session_state.gene_exp)

    # Find compounds
    st.subheader("Route search")
    cpd_query = st.text_input("Enter compound keyword:")
    if cpd_query:
        biocyc_ids = cpds.search_compound_biocyc_id(cpd_query)
        if debug:
            st.write(biocyc_ids)

        cpd_id = st.radio("Compound ID to use:", biocyc_ids["hits"], horizontal=True)
        compartments = cpds.search_compound_compartment(cpd_id)
        compartment = st.radio("Cellular compartment:", compartments, horizontal=True)
        cpd_metaid = cpds.search_compound_metaid_in_compartment(cpd_id, compartment)

        if st.button("Route search"):
            with st.spinner(f"Searching routes starting from {cpd_metaid}"):
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
                num_routes = min(len(routes), max_num_routes)
                st.success(
                    f"Found {len(routes)} routes in total, returning top {num_routes}."
                )
                for r in routes:
                    r["genes"] = "|".join(list(r["genes"]))

                res = (
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
                    .head(max_num_routes)
                )

                st.dataframe(res, use_container_width=True)

                # Provide download link
                csv = download_df(res)
                st.download_button(
                    label="Download data as CSV",
                    data=csv,
                    file_name=f"{cpd_id}-route-search.csv",
                    mime="text/csv",
                )

con.close()
