# Overview

![metabolike logo](_static/metabolike-logo-round.png){ align=right }

**An alternative way to explore metabolic networks.**

Metabolike lets you transform SBML metabolic models into queryable, interactive graphs.
It provides a unified interface towards studying metabolic reprogramming under perturbed cellular states.

At its core, metabolike consists of two parts:
a parser that can collect data from SBML models and transform the information into a graph,
and a set of graph algorithms that powers the reprogrammed route search within the graph.
Metabolike also includes a [Streamlit][api.main] app for easy route search with both example TCGA datasets and custom CSV files.

To get started, follow [installation instructions](usage/installation.md) to get the CLI tool.
Then you can [start building your first metabolic graph](usage/getting_started.md)!
