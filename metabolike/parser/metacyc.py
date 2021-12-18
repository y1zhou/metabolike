import logging
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
]


class Metacyc:
    """
    Converting MetaCyc files to a Neo4j database.
    Documentation on the MetaCyc files and ofrmat FAQs can be found at:

    - MetaCyc data files download: https://metacyc.org/downloads.shtml
    - MetaCyc file formats: http://bioinformatics.ai.sri.com/ptools/flatfile-format.html
    - SBML FAQ: https://synonym.caltech.edu/documents/faq
    """

    def __init__(self, filepath: Union[str, Path], neo4j_driver: db.Neo4jDriver):
        self.neo4j_driver = neo4j_driver
        self.filepath = Path(filepath).expanduser().resolve()

    def setup(self, force: bool = False):
        # Read SBML file
        self.doc = self.read_sbml()
        self.model: libsbml.Model = self.doc.getModel()
        self.db_name: str = self.model.getMetaId().lower()

        # Setup Neo4j database
        self.setup_graph_db(force=force)

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

    def setup_graph_db(self, **kwargs):
        """
        Create Neo4j database and set proper constraints. `Reaction` nodes are
        central in the schema:

        ```
                     Pathway
                        │
        Compartment───Reaction───Compound
            │           │          │
            └───────────┼──────────┘
                        │
                    GeneProduct
        ```
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
