# Overview

![metabolike logo](_static/metabolike-logo-round.png){ align=right }

This is a Python package that aggregates the data from
[BRENDA](https://brenda-enzymes.org/) and [BioCyc](https://biocyc.org/)
into one unified [graph database](https://neo4j.com/). The MetaCyc
(BioCyc) database is imported to a graph database using the provided
SBML file and various `.dat` files. A [Streamlit](https://streamlit.io/)
app is [included in the
code](https://github.com/y1zhou/metabolike/blob/main/metabolike/api/main.py)
for easy route search with both example TCGA datasets and custom CSV
files.
