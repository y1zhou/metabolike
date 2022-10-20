Metabolic reprogramming Knowledgebase
=====================================

.. image:: https://img.shields.io/pypi/v/metabolike
    :alt: PyPI release
    :target: https://pypi.org/project/metabolike/
.. image:: https://readthedocs.org/projects/metabolike/badge/?version=stable&style=flat
    :alt: Documentation status for stable version
    :target: https://metabolike.readthedocs.io/en/stable/
.. image:: https://img.shields.io/github/workflow/status/y1zhou/metabolike/PyTest?label=test
    :alt: PyTest workflow status
    :target: https://github.com/y1zhou/metabolike/actions/workflows/pytest.yml
.. image:: https://codecov.io/gh/y1zhou/metabolike/branch/main/graph/badge.svg
    :alt: Code coverage
    :target: https://codecov.io/gh/y1zhou/metabolike
.. image:: https://github.com/y1zhou/metabolike/workflows/CodeQL/badge.svg
    :alt: CodeQL analysis
    :target: https://github.com/y1zhou/metabolike/actions/workflows/codeql-analysis.yml
.. image:: https://img.shields.io/github/license/y1zhou/metabolike
    :alt: Package license
    :target: https://github.com/y1zhou/metabolike/blob/main/LICENSE
.. image:: https://img.shields.io/badge/code%20style-black-000000
    :alt: Code styled with Black
    :target: https://github.com/psf/black
.. image:: https://img.shields.io/github/commits-since/y1zhou/metabolike/latest/main
    :alt: GitHub commits since latest release (by date) for main branch
    :target: https://github.com/y1zhou/metabolike/commits/main

A Python package that aggregates the data from BRENDA_ and BioCyc_ into one unified `graph database`_.
The MetaCyc (BioCyc) database is imported to a graph database using the provided SBML file and various ``.dat`` files.
A `Streamlit`_ app is `included in the code`_ for easy route search with both example TCGA datasets and custom CSV files.

.. _BRENDA: https://brenda-enzymes.org/
.. _BioCyc: https://biocyc.org/
.. _graph database: https://neo4j.com/
.. _Streamlit: https://streamlit.io/
.. _included in the code: https://github.com/y1zhou/metabolike/blob/main/metabolike/api/main.py

To-Do
-----

* Add tutorial for the package
* Add Dockerfile and docker-compose file
