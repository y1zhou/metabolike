import os

import pandas as pd
import streamlit as st
from streamlit_agraph import Edge, Node

from metabolike.algo.compounds import CompoundMap
from metabolike.algo.omics import ReactionGeneMap
from metabolike.db import Neo4jClient

TCGA = (
    "TCGA-BLCA",
    "TCGA-BRCA",
    "TCGA-COAD",
    "TCGA-ESCA",
    "TCGA-HNSC",
    "TCGA-KICH",
    "TCGA-KIRC",
    "TCGA-KIRP",
    "TCGA-LIHC",
    "TCGA-LUAD",
    "TCGA-LUSC",
    "TCGA-PRAD",
    "TCGA-STAD",
    "TCGA-THCA",
)


@st.experimental_singleton
def neo4j_db():
    """Connects to Neo4j database."""
    neo4j_uri = os.environ.get("NEO4J_URI", "neo4j://localhost:7687")
    neo4j_password = os.environ.get("NEO4J_PASSWORD")
    neo4j_database = os.environ.get("NEO4J_DATABASE", "neo4j")
    _con = Neo4jClient(
        uri=neo4j_uri,
        neo4j_password=neo4j_password,
        database=neo4j_database,
    )
    return _con


@st.experimental_singleton
def get_compound_map(_con: Neo4jClient):
    return CompoundMap(_con)


@st.experimental_singleton
def download_tcga_data(tcga_project: str):
    tcga_url = f"https://s3.y1zhou.com/{tcga_project}.csv"
    req_headers = {
        "User-Agent": "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.11 (KHTML, like Gecko) Chrome/23.0.1271.64 Safari/537.11",
        "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
        "Accept-Charset": "ISO-8859-1,utf-8;q=0.7,*;q=0.3",
        "Accept-Encoding": "none",
        "Accept-Language": "en-US,en;q=0.8",
        "Connection": "keep-alive",
    }
    df = pd.read_csv(tcga_url, storage_options=req_headers)
    return df


@st.cache
def download_df(df):
    # IMPORTANT: Cache the conversion to prevent computation on every rerun
    return df.to_csv().encode("utf-8")


def routes_to_graph(
    routes, rgm: ReactionGeneMap, pos_color="#FF4B4B", neg_color="#002941"
):
    """
    Generates a graph object used by streamlit for visualizing the routes.
    Duplicated reaction nodes are dropped but edges are preserved.
    Edge widths are proportional to the highest-scoring route they belong to.
    """
    nodes, edges = {}, {}
    max_score, min_score = -float("inf"), float("inf")

    for route_id, route in enumerate(routes):
        # Save max/min scores for later scaling
        max_score = max(max_score, abs(route["score"]))
        min_score = min(min_score, abs(route["score"]))

        # Iterate over all nodes
        r = route["route"]
        route_color = pos_color if route["score"] > 0 else neg_color
        for i, step in enumerate(r):
            # compound->reaction->compound|compounds->reaction->...
            if step["nodeType"] == "Compound":
                # Default state for new nodes
                if step["id"] not in nodes:
                    n = Node(
                        id=step["id"],
                        label=step["name"],
                        size=6,
                        shape="dot",
                        color="#EA9C70",
                    )
                    nodes[step["id"]] = n
                # TODO: color nodes based on structural similarity?

            elif step["nodeType"] == "Reaction":
                n_src = r[i - 1]["id"]
                _ = edges.setdefault(n_src, {})

                k = i + 1
                while k < len(r) and r[k]["nodeType"] == "Compound":
                    n_tgt = r[k]["id"]

                    if n_tgt not in edges[n_src]:
                        rxn_exp = rgm.rxn_exp.get(step["id"])
                        if rxn_exp is None:
                            rxn_exp = route["score"]
                        edges[n_src][n_tgt] = Edge(
                            source=n_src,
                            target=n_tgt,
                            title=f"route {route_id}: {step['name']}",
                            value=abs(rxn_exp),
                            color=route_color,
                        )

                    k += 1

            else:
                raise ValueError(f"Unknown node type {step['nodeType']}")

    # Scale edge widths
    final_edges = [x for n in edges.values() for x in n.values()]
    scaling_factor = max_score - min_score
    for e in final_edges:
        e.value = (e.value - min_score) / scaling_factor
        e.value = min(e.value, 0.75)

    return list(nodes.values()), final_edges
