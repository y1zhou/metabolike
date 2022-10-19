import os

import pandas as pd
import streamlit as st

from metabolike.algo.compounds import CompoundMap
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
