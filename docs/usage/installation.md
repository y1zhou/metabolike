# Installation

## Prerequisites

You will need a running Neo4j database and login credentials of a user with write access to the database.
To populate the database, you will also need to acquire the [BioCyc data files](https://biocyc.org/download.shtml) and the [BRENDA text file](https://www.brenda-enzymes.org/download_brenda_without_registration.php).
The `metabolike` package doesn't include the data sources since they are large and require license agreements. Please download the data files and extract them for later use.

The MetaCyc (BioCyc) database is imported to the graph database using the provided SBML file and various `.dat` files.
The BRENDA text file is then mapped onto the graph using common EC numbers and KEGG reaction IDs.

## Latest stable version

To install the latest stable version of the package, run the following command:

```bash
pip install metabolike
```

## Development version

As the package is still in its early stages, it is recommended to get the latest version from GitHub:

```bash
git clone --depth 1 https://github.com/y1zhou/metabolike
```

To install the package, run the following command:

```bash
pip install . # (1)
```

1. Note that you should `cd` into the directory first, i.e. `#!bash cd metabolike`.

You can also install the package using [flit](https://flit.readthedocs.io/en/latest/):

```shell
flit install
```
