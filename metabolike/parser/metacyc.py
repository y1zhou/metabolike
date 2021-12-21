import logging
from itertools import groupby
from pathlib import Path
from typing import Dict, List, Tuple, Union

import libsbml
from metabolike import db

logger = logging.getLogger(__name__)

# MIRIAM qualifiers (https://co.mbine.org/standards/qualifiers)
# libsbml.BiolQualifierType_toString
BIO_QUALIFIERS = {
    0: "is",
    1: "hasPart",
    2: "isPartOf",
    3: "isVersionOf",
    4: "hasVersion",
    5: "isHomologTo",
    6: "isDescribedBy",
    7: "isEncodedBy",
    8: "encodes",
    9: "occursIn",
    10: "hasProperty",
    11: "isPropertyOf",
    12: "hasTaxon",
    13: "unknown",
}

NODE_LABELS = [
    "Pathway",
    "Compartment",
    "Reaction",
    "Compound",
    "GeneProduct",
    "RDF",
    "Complex",
    "EntitySet",
]


class Metacyc:
    """
    Converting MetaCyc files to a Neo4j database.
    Documentation on the MetaCyc files and ofrmat FAQs can be found at:

    - MetaCyc data files download: https://metacyc.org/downloads.shtml
    - MetaCyc file formats: http://bioinformatics.ai.sri.com/ptools/flatfile-format.html
    - SBML FAQ: https://synonym.caltech.edu/documents/faq
    """

    def __init__(
        self, filepath: Union[str, Path], neo4j_driver: db.Neo4jDriver, **kwargs
    ):
        """
        Args:
            filepath: The path to the MetaCyc SBML file to convert.
            neo4j_driver: A Neo4jDriver instance.
            pw_filepath (optional): The path to the pathway.dat file. If given,
                the file will be parsed and pathway links will be added to the
                Reaction nodes.
        """
        # File paths
        self.filepath = self._validate_path(filepath)
        pw_filepath = kwargs.get("pw_filepath")
        if pw_filepath:
            self.pw_filepath = self._validate_path(pw_filepath)

        # Neo4j driver
        self.neo4j_driver = neo4j_driver

    def setup(self, force: bool = False):
        # Read SBML file
        self.doc = self.read_sbml()
        self.model: libsbml.Model = self.doc.getModel()
        self.db_name: str = self.model.getMetaId().lower()

        # Setup Neo4j database and populate it with data
        self.setup_graph_db(force=force)
        self.sbml_to_graph()

        # Read pathways file if given
        if self.pw_filepath:
            # `rb` to bypass the 0xa9 character
            with open(self.pw_filepath, "rb") as f:
                lines = [l.decode("utf-8", "ignore").strip() for l in f]
                # Also remove comments on the top of the file
                lines = [l for l in lines if not l.startswith("#")]

            # Split entries based on `//`
            doc = [list(g) for k, g in groupby(lines, key=lambda x: x != "//") if k]
            doc = [x for x in doc if len(x) > 1]  # Remove empty entries
            pathways = {}
            for pw in doc:
                pw_id, pw_data = self.read_pathways(pw)
                pathways[pw_id] = pw_data

    def read_sbml(self) -> libsbml.SBMLDocument:
        reader = libsbml.SBMLReader()
        metacyc = reader.readSBMLFromFile(self.filepath)
        logger.info("Finished reading SBML file")

        for i in range(metacyc.getNumErrors()):
            err = metacyc.getError(i)
            logger.warning(
                "SBML reader raised %s during parsing at line %d, column %d: %s",
                err.getSeverityAsString(),
                err.getLine(),
                err.getColumn(),
                err.getMessage().replace("\n", " "),
            )

        return metacyc

    @staticmethod
    def read_pathways(lines: List[str]) -> Tuple[str, Dict[str, Union[str, List[str]]]]:
        """
        Parse one entry from an attribute-value file.

        Args:
            lines: A list of lines from an entry of the file.
        """
        pw_id, pathway = "", {}
        for line in lines:
            # / marks the continuation of the previous line
            if line.startswith("/"):
                # TODO
                continue

            # ^ marks an annotation-value pair where the annotation is a
            # label attached to the previous attribute
            elif line.startswith("^"):
                # TODO
                continue

            # Otherwise we have a regular attribute-value pair
            else:
                attr, val = line.split(" - ")
                # Treat `UNIQUE-ID` as a special case as we need to return
                # the pathway ID separately
                if attr == "UNIQUE-ID":
                    pw_id = val
                    continue
                # TODO

        if not pw_id:
            raise ValueError("Pathway ID not found")

        return pw_id, pathway

    def setup_graph_db(self, **kwargs):
        """
        Create Neo4j database and set proper constraints. `Reaction` nodes are
        central in the schema (see README.rst).
        """
        # Create database
        try:
            db.create(self.neo4j_driver, self.db_name, **kwargs)
        except Exception as e:
            logger.fatal(f"Could not create database: {e}")
            raise

        # Set constraints
        with self.neo4j_driver.session(database=self.db_name) as session:
            logger.debug("Creating constraint for RDF nodes")
            r = session.run(
                """CREATE CONSTRAINT IF NOT EXISTS
                   ON (r:RDF) ASSERT r.uri IS UNIQUE;"""
            ).data()
            if r:
                logger.warning(f"Could not create constraint for RDF nodes: {r}")

            # Constraints automatically create indexes, so we don't need to
            # create them manually.
            for label in NODE_LABELS:
                logger.debug(f"Creating constraint for {label} nodes")
                r = session.run(
                    f"""CREATE CONSTRAINT IF NOT EXISTS
                        ON (n:{label}) ASSERT n.mcId IS UNIQUE;""",
                ).data()
                if r:
                    logger.warning(
                        f"Could not create constraint for {label} nodes: {r}"
                    )

    def sbml_to_graph(self):
        """
        Populate Neo4j database with SBML data.

        Nodes are created for each SBML element using `MERGE` statements:
        https://neo4j.com/docs/cypher-manual/current/clauses/merge/#merge-merge-with-on-create
        """
        with self.neo4j_driver.session(database=self.db_name) as session:
            # Compartments
            logger.info("Creating Compartment nodes")
            for c in self.model.getListOfCompartments():
                c: libsbml.Compartment
                session.write_transaction(
                    lambda tx: tx.run(
                        """
                        MERGE (c:Compartment {mcId: $mcId})
                        ON CREATE
                        SET c.displayName = $name;
                        """,
                        mcId=c.getMetaId(),
                        name=c.getName(),
                    ),
                )

            # Compounds, i.e. metabolites, species
            logger.info("Creating Compound nodes")
            for s in self.model.getListOfSpecies():
                s: libsbml.Species
                # Basic properties
                mcid: str = s.getMetaId()
                logger.debug(f"Creating Compound node {mcid}")
                session.write_transaction(
                    lambda tx: tx.run(
                        """
                        MATCH (cpt:Compartment {mcId: $compartment})
                        MERGE (c:Compound {mcId: $mcId})-[:hasCompartment]->(cpt)
                        ON CREATE
                            SET c.displayName = $name,
                                c.charge = $charge,
                                c.chemicalFormula = $formula,
                                c.boundaryCondition = $boundaryCondition,
                                c.hasOnlySubstanceUnits = $hasOnlySubstanceUnits,
                                c.constant = $constant;
                        """,
                        mcId=mcid,
                        name=s.getName(),
                        compartment=s.getCompartment(),
                        charge=s.getCharge(),
                        formula=s.getPlugin("fbc").getChemicalFormula(),
                        boundaryCondition=s.getBoundaryCondition(),  # unchanged by reactions
                        hasOnlySubstanceUnits=s.getHasOnlySubstanceUnits(),
                        constant=s.getConstant(),
                    )
                )

                # Add RDF annotations
                logger.debug(f"Adding RDF nodes for Compound {mcid}")
                cvterms: List[libsbml.CVTerm] = s.getCVTerms()
                for cvterm in cvterms:
                    self.add_rdf_node("Compound", mcid, cvterm, session)

            # Gene products
            logger.info("Creating GeneProduct nodes")
            fbc: libsbml.FbcModelPlugin = self.model.getPlugin("fbc")
            for gp in fbc.getListOfGeneProducts():
                gp: libsbml.GeneProduct
                mcid = gp.getMetaId()
                logger.debug(f"Creating GeneProduct node {mcid}")
                session.write_transaction(
                    lambda tx: tx.run(
                        """
                        MERGE (gp:GeneProduct {mcId: $mcId})
                        ON CREATE
                            SET gp.displayName = $name,
                                gp.label = $label;
                        """,
                        mcId=mcid,
                        name=gp.getName(),
                        label=gp.getLabel(),
                    ),
                )
                logger.debug(f"Adding RDF nodes for GeneProduct {mcid}")
                cvterms: List[libsbml.CVTerm] = gp.getCVTerms()
                for cvterm in cvterms:
                    self.add_rdf_node("GeneProduct", mcid, cvterm, session)

            # Reactions
            logger.info("Creating Reaction nodes")
            for r in self.model.getListOfReactions():
                r: libsbml.Reaction
                mcid = r.getMetaId()

                # Basic properties
                logger.debug(f"Creating Reaction node {mcid}")
                session.write_transaction(
                    lambda tx: tx.run(
                        """
                        MERGE (r:Reaction {mcId: $mcId})
                        ON CREATE
                            SET r.displayName = $name,
                                r.reversible = $reversible,
                                r.fast = $fast;
                        """,
                        mcId=mcid,
                        name=r.getName(),
                        reversible=r.getReversible(),
                        fast=r.getFast(),
                    ),
                )

                # Add RDF nodes
                logger.debug(f"Adding RDF nodes for Reaction {mcid}")
                cvterms: List[libsbml.CVTerm] = r.getCVTerms()
                for cvterm in cvterms:
                    self.add_rdf_node("Reaction", mcid, cvterm, session)

                # Add reactants and products
                logger.debug(f"Adding reactants for Reaction {mcid}")
                reactants = r.getListOfReactants()
                self.link_reaction_to_compound(mcid, reactants, "Reactant", session)

                logger.debug(f"Adding products for Reaction {mcid}")
                products = r.getListOfProducts()
                self.link_reaction_to_compound(mcid, products, "Product", session)

                # Add associated gene products
                # This could be complicated where the child nodes could be:
                # 1. GeneProductRef
                # 2. fbc:or -> GeneProductRef
                # 3. fbc:and -> GeneProductRef
                # 4. Arbitrarily nested fbc:or / fbc:and within cases 2 and 3
                gpa: libsbml.GeneProductAssociation = r.getPlugin(
                    "fbc"
                ).getGeneProductAssociation()
                if gpa is not None:
                    node = gpa.getAssociation()
                self.add_gene_product_association_node(node, session, mcid)

    def add_rdf_node(
        self, node_label: str, mcid: str, cvterm: libsbml.CVTerm, session: db.Session
    ):
        """Create RDF node and link it to the given SBML node.

        RDF in Annotation are in the form of triples:
        the model component to annotate (subject), the relationship between
        the model component and the annotation (predicate), and a term
        describing the component (object).

        Args:
            node_label: The label of the SBML node to annotate. Should be one
                of the labels defined in `NODE_LABELS`.
            mcid: The MetaCyc ID of the SBML node to annotate.
            cvterm: The CVTerm that contains information to add to the RDF node.
            session: The Neo4j session to use.
        """
        # Get the biological qualifier type of the terms
        bio_qual = BIO_QUALIFIERS[cvterm.getBiologicalQualifierType()]
        # Get the content of each RDF term
        uris = [
            self.split_uri(cvterm.getResourceURI(i))
            for i in range(cvterm.getNumResources())
        ]
        uris = {x[0]: x[1] for x in uris}
        # Generate the Cypher statements for each RDF term
        uri_cypher = self.generate_cypher_for_uri(uris)

        session.write_transaction(
            lambda tx: tx.run(
                f"""
                MATCH (c:{node_label} {{mcId: $mcId}})
                MERGE (n:RDF)<-[:{bio_qual}]-(c)
                ON CREATE
                    SET {uri_cypher};
                """,
                mcId=mcid,
                **uris,
            )
        )

    @staticmethod
    def split_uri(uri: str) -> Tuple[str, str]:
        """Split a URI into a namespace and an annotation term.

        Args:
            uri: URI to split.

        Returns:
            Tuple of namespace and annotation term.
        """
        # First three elements are from http://identifiers.org/
        uri = uri.replace("-", ".")  # Ec-code
        res = uri.split("/")[3:]
        resource, identifier = res[0], res[1:]
        resource = "".join(x.capitalize() for x in resource.split("."))
        identifier = "".join(identifier)
        return resource, identifier

    @staticmethod
    def generate_cypher_for_uri(uris: Dict[str, str]) -> str:
        cypher = []
        for label in uris.keys():
            cypher.append(f"n.{label} = ${label}")

        return ",".join(cypher)

    @staticmethod
    def link_reaction_to_compound(
        reaction_id: str,
        compounds: List[libsbml.SpeciesReference],
        compound_type: str,
        session: db.Session,
    ):
        """Link reactants or products to a reaction.

        Args:
            reaction_id: The MetaCyc ID of the reaction.
            compounds: The list of compounds to link to the reaction.
            compound_type: The type of compound to link to the reaction. Should
                be one of "Reactant" or "Product".
            session: The Neo4j session to use.
        """
        if compound_type not in ["Reactant", "Product"]:
            raise ValueError(f"Invalid compound type: {compound_type}")
        for cpd in compounds:
            logger.debug(
                f"Adding {compound_type} {cpd.getSpecies()} to Reaction {reaction_id}"
            )
            session.write_transaction(
                lambda tx: tx.run(
                    f"""
                    MATCH (r:Reaction {{mcId: $reaction}}),
                            (c:Compound {{mcId: $compound}})
                    MERGE (r)-[l:has{compound_type}]->(c)
                    ON CREATE
                        SET l.stoichiometry = $stoichiometry,
                            l.constant = $constant;
                    """,
                    reaction=reaction_id,
                    compound=cpd.getSpecies(),
                    stoichiometry=cpd.getStoichiometry(),
                    constant=cpd.getConstant(),
                )
            )

    def add_gene_product_association_node(
        self,
        node: Union[libsbml.GeneProductRef, libsbml.FbcAnd, libsbml.FbcOr],
        session: db.Session,
        source_id: str,
        source_label: str = "Reaction",
        edge_type: str = "hasGeneProduct",
        node_index: int = 0,
    ):
        """
        Add gene products to a reaction. When the added node is FbcAnd or FbcOr,
        recursively add the children. This means a custom `mcId` is constructed
        for the `Complex` and `EntitySet` nodes corresponding to the `FbcAnd`
        and `FbcOr` nodes, respectively.

        Args:
            node: The GeneProductAssociation child node to add.
            session: The Neo4j session to use.
            source_id: The MetaCyc ID of the source node. This should be the
                MetaCyc ID of the `Reaction` node
            source_label: The label of the source node.
            edge_type: The type of edge to add. Should be one of `hasGeneProduct`,
                `hasComponent`, or `hasMember`.
            node_index: The index of the current node. This is used to construct
                the `mcId` of the `Complex` and `EntitySet` nodes.
        """
        # If there's no nested association, add the node directly
        if isinstance(node, libsbml.GeneProductRef):
            session.write_transaction(
                lambda tx: tx.run(
                    f"""
                    MATCH (n:{source_label} {{mcId: $node_id}}),
                          (gp:GeneProduct {{mcId: $gp_id}})
                    MERGE (gp)<-[l:{edge_type}]-(n)
                    """,
                    node_id=source_id,
                    gp_id=node.getGeneProduct(),
                )
            )

        # For nested associations, first add a `Complex` or `EntitySet` node,
        # then recursively add the children
        elif isinstance(node, libsbml.FbcAnd):
            complex_id = f"{source_id}_complex{node_index}"
            session.write_transaction(
                lambda tx: tx.run(
                    f"""
                    MATCH (n:{source_label} {{mcId: $node_id}})
                    MERGE (:Complex {{mcId: $complex_id}})<-[l:{edge_type}]-(n)
                    """,
                    node_id=source_id,
                    complex_id=complex_id,
                )
            )
            for i in range(node.getNumAssociations()):
                self.add_gene_product_association_node(
                    node.getAssociation(i),
                    session,
                    complex_id,
                    "Complex",
                    "hasComponent",
                    i,
                )
        elif isinstance(node, libsbml.FbcOr):
            eset_id = f"{source_id}_entityset{node_index}"
            session.write_transaction(
                lambda tx: tx.run(
                    f"""
                    MATCH (n:{source_label} {{mcId: $node_id}})
                    MERGE (:EntitySet {{mcId: $eset_id}})<-[l:{edge_type}]-(n)
                    """,
                    node_id=source_id,
                    eset_id=eset_id,
                )
            )
            for i in range(node.getNumAssociations()):
                self.add_gene_product_association_node(
                    node.getAssociation(i),
                    session,
                    eset_id,
                    "EntitySet",
                    "hasMember",
                    i,
                )
        else:
            logging.error(f"Unhandled GeneProductAssociation type {type(node)}")
            raise ValueError

    @staticmethod
    def _validate_path(filepath: Union[str, Path]) -> Path:
        f = Path(filepath).expanduser().resolve()
        if not f.is_file():
            logger.error(f"File does not exist: {f}")
            raise FileNotFoundError(str(f))

        return f
