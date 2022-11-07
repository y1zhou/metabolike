[![PyPI release](https://img.shields.io/pypi/v/metabolike)](https://pypi.org/project/metabolike/)
[![Documentation status for stable version](https://readthedocs.org/projects/metabolike/badge/?version=stable&style=flat)](https://metabolike.readthedocs.io/en/stable/)
[![PyTest workflow status](https://img.shields.io/github/workflow/status/y1zhou/metabolike/PyTest?label=test)](https://github.com/y1zhou/metabolike/actions/workflows/pytest.yml)
[![Code coverage](https://codecov.io/gh/y1zhou/metabolike/branch/main/graph/badge.svg)](https://codecov.io/gh/y1zhou/metabolike)
[![CodeQL analysis](https://github.com/y1zhou/metabolike/workflows/CodeQL/badge.svg)](https://github.com/y1zhou/metabolike/actions/workflows/codeql-analysis.yml)
[![Package license](https://img.shields.io/github/license/y1zhou/metabolike)](https://github.com/y1zhou/metabolike/blob/main/LICENSE)
[![Code styled with Black](https://img.shields.io/badge/code%20style-black-000000)](https://github.com/psf/black)
[![GitHub commits since latest release (by date) for main branch](https://img.shields.io/github/commits-since/y1zhou/metabolike/latest/main)](https://github.com/y1zhou/metabolike/commits/main)

This is a Python package that aggregates the data from
[BRENDA](https://brenda-enzymes.org/) and [BioCyc](https://biocyc.org/)
into one unified [graph database](https://neo4j.com/). The MetaCyc
(BioCyc) database is imported to a graph database using the provided
SBML file and various `.dat` files. A [Streamlit](https://streamlit.io/)
app is [included in the
code](https://github.com/y1zhou/metabolike/blob/main/metabolike/api/main.py)
for easy route search with both example TCGA datasets and custom CSV
files.

For more details, please visit the [documentation site](https://metabolike.readthedocs.io/en/stable/).
