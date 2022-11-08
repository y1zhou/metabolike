# Getting started

## Prerequisites

In this section, you will make use of the parser module of `metabolike`.
More specifically, you will use the \[`SBML parser module`\]\[parser.sbml\] to extract information from an SBML model file and populate the graph database.
But before that, you will need to set up Neo4j and download a sample SBML file.

### Neo4j database

Since the main function of this package is to convert metabolic networks from text files to graphs, a running [Neo4j database](https://neo4j.com/product/neo4j-graph-database/) instance is needed as your storage backend.
You also need login credentials of a user with **write access** to the database.

If you prefer a managed cloud database, the official [Neo4j AuraDB](https://neo4j.com/cloud/platform/aura-graph-database/) is a good resource.
Once you set up an account with them, follow the onboarding process and create an **empty** database instance.
You will be prompted with a username and default password.
After a few minutes, the new instance will be created and you will get a connection URI in the form of `neo4j+s://<id>.databases.neo4j.io`.
Keep track of these three values as you will need them soon!
You may now jump to the [SBML model files](#sbml-model-files) section.

#### Using Docker

!!! tip

    This assumes Docker is already up and running on your machine.
    If not, no worries! [Installation of Docker](https://docs.docker.com/engine/install/) is a one-liner on most systems.
    More specific configurations could take a _lot_ longer, but I digress.

The other way to go is to set up your own Neo4j database.
Installation of Neo4j is explained in great detail on the [official documentation website](https://neo4j.com/docs/operations-manual/current/installation/).
For simplicity, this tutorial will focus on the Docker setup.

The Neo4j team maintains [an official image on DockerHub](https://hub.docker.com/_/neo4j) and has [great docs about its usage](https://neo4j.com/docs/operations-manual/current/docker/introduction/).
The following is sufficient to spin up a Neo4j instance listening on [localhost:7474](http://localhost:7474):

```bash
docker run \
    --restart always \
    -p=7474:7474 \
    -p=7687:7687 \ # (1)
    neo4j:4.4.9 # (2)
```

1. Port `7687` is used by the Bolt protocol.
   Check [the manual](https://neo4j.com/docs/operations-manual/current/configuration/ports/) or all ports relevant to Neo4j.
2. If you have an enterprise license, use `neo4j:4.4.9-enterprise` instead.

You most likely want to persist data across restarts, so Docker volumes should also be added:

```bash hl_lines="5 6"
docker run \
    --restart always \
    -p=7474:7474 \
    -p=7687:7687 \
    -v /path/to/data:/data \
    -v /path/to/logs:/logs \ # (1)
    neo4j:4.4.9
```

1. The mounted directories `/path/to/[data|logs]` should exist **before** running Docker, otherwise you'll either run
   into permission errors or have files created instead of directories.

The `metabolike` package requires the [APOC library](https://neo4j.com/labs/apoc/4.3/installation/) for complex queries.
You will have to download the binary `jar` file from [their release page](http://github.com/neo4j-contrib/neo4j-apoc-procedures/releases/) and place it in a directory (say `plugins/`).
For Neo4j 4.4.9, the APOC file you need is `apoc-4.3.0.9-all.jar`.
When this directory is mounted, the Neo4j image automatically recognizes and loads the plugin:

```bash hl_lines="7"
docker run \
    --restart always \
    -p=7474:7474 \
    -p=7687:7687 \
    -v /path/to/data:/data \
    -v /path/to/logs:/logs \
    -v /path/to/plugins:/plugins \
    neo4j:4.4.9
```

#### Using Docker Compose

If you happen to have [Docker Compose](https://docs.docker.com/compose/) installed, grab the following file:

```yaml title="docker-compose.yaml"
version: "3"

services:
  neo4j:
    image: "neo4j:4.4.9"
    container_name: neo4j-metabolike
    restart: unless-stopped
    ports:
      - 7474:7474
      - 7687:7687
    volumes:
      - neo4j-data:/data
      - neo4j-logs:/logs
      - "${PWD}/plugins:/plugins" # (1)
    environment: # (2)
      - NEO4J_apoc_export_file_enabled=true
      - NEO4J_apoc_import_file_enabled=true
      - NEO4J_apoc_import_file_use__neo4j__config=false
      - NEO4J_dbms_allow__upgrade=true
      - NEO4J_dbms_security_procedures_unrestricted=gds.*,apoc.*
      - NEO4J_dbms_security_procedures_allowlist=gds.*,apoc.*
      - NEO4J_ACCEPT_LICENSE_AGREEMENT=yes # (3)
    env_file:
      - .env # (4)

volumes:
  neo4j-data:
  neo4j-logs:
```

1. You should have a `plugins/` directory alongside `docker-compose.yaml`, and `apoc-4.3.0.9-all.jar` in the directory.
2. The environment variables are required by APOC in some cases.
3. Include this line if you are using the enterprise licensed version.
4. You can pass additional environment variables in a `.env` file in the same directory.

Inside the directory with the YAML file, run:

```bash
docker compose up -d
```

Now you can head to [localhost:7474](http://localhost:7474) and login using the default username and password.
You may be prompted to change the password after the first successful connection.

??? info "Default credentials"

    By default, the initial values for both the username and password are `neo4j`.
    You can override this by providing an environment variable `NEO4J_AUTH=neo4j/your-password`.
    You can also set `NEO4J_AUTH=none` and disable authentication for local debugging purposes.

### SBML model files

### Optional files

#### MetaCyc

#### BRENDA

To populate the database, you will also need to acquire the [BioCyc data files](https://biocyc.org/download.shtml) and
the [BRENDA text file](https://www.brenda-enzymes.org/download_brenda_without_registration.php).
The `metabolike` package doesn't include the data sources since they are large and require license agreements. Please
download the data files and extract them for later use.

The MetaCyc (BioCyc) database is imported to the graph database using the provided SBML file and various `.dat` files.
The BRENDA text file is then mapped onto the graph using common EC numbers and KEGG reaction IDs.

## SBML model parsing

A schema of the final built graph is:

![Schema of the MetaCyc database](../_static/metabolike_schema.svg)

<figcaption>Schema of the MetaCyc database</figcaption>
