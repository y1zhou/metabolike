# Welcome to metabolike :wave:

<img align="right" width="150" height="150" src="./docs/_static/metabolike-logo-round.png" alt="Metabolike logo">

[![PyPI release](https://img.shields.io/pypi/v/metabolike)](https://pypi.org/project/metabolike/)
[![Documentation status for stable version](https://readthedocs.org/projects/metabolike/badge/?version=stable&style=flat)](https://metabolike.readthedocs.io/en/stable/)
[![PyTest workflow status](https://img.shields.io/github/workflow/status/y1zhou/metabolike/PyTest?label=test)](https://github.com/y1zhou/metabolike/actions/workflows/pytest.yml)
[![Code coverage](https://codecov.io/gh/y1zhou/metabolike/branch/main/graph/badge.svg)](https://codecov.io/gh/y1zhou/metabolike)
[![CodeQL analysis](https://github.com/y1zhou/metabolike/workflows/CodeQL/badge.svg)](https://github.com/y1zhou/metabolike/actions/workflows/codeql-analysis.yml)
[![Package license](https://img.shields.io/github/license/y1zhou/metabolike)](https://github.com/y1zhou/metabolike/blob/main/LICENSE)
[![Code styled with Black](https://img.shields.io/badge/code%20style-black-000000)](https://github.com/psf/black)
[![GitHub commits since latest release (by date) for main branch](https://img.shields.io/github/commits-since/y1zhou/metabolike/latest/main)](https://github.com/y1zhou/metabolike/commits/main)

**An alternative way to explore metabolic networks.**

Metabolike lets you transform SBML metabolic models into queryable, interactive graphs.

## Installation

```bash
pip install -U metabolike
```

For more details, please visit the [documentation site](https://metabolike.rtfd.io/).

## Feature highlights

- Parser that can read any valid SBML model and write relevant entities into a [graph database](https://neo4j.com/).
- Special support for [BRENDA](https://brenda-enzymes.org/) and [BioCyc](https://biocyc.org/) data files.
- A [Streamlit](https://streamlit.io/)
  app for novel metabolic route detection with omics-data integration.

## License

Metabolike is completely free and open-source and licensed under the GPL-3.0 license.
