# Getting started

The project layout is roughly as follows:

```yaml
metabolike: # (1)
  algo:
    - compounds.py
    - omics.py
    - routes.py
  api:
    - data.py
    - main.py # (2)
  db:
    - metacyc.py
    - neo4j.py
    - sbml.py
  parser:
     - bkms_react.py
     - brenda.py
     - brenda_transformer.py
     - metacyc.py
     - sbml.py
```

1. Top-level models include `config.py`, `main.py`, and `utils.py`.
2. The streamlit app entrypoint.

A schema of the final built graph is:

![Schema of the MetaCyc database](../_static/metabolike_schema.svg)

<figcaption>Schema of the MetaCyc database</figcaption>
