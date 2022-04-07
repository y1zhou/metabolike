import logging
import re
from pathlib import Path
from typing import Dict, List, Tuple, Union

import libsbml
from metabolike.db import SBMLClient
from metabolike.utils import add_kv_to_dict, validate_path
from tqdm import tqdm

logger = logging.getLogger(__name__)


class SBMLParser:
    """
    Converting MetaCyc files to a Neo4j database.
    Documentation on the MetaCyc files and format FAQs can be found at:

    * MetaCyc data files download: https://metacyc.org/downloads.shtml
    * MetaCyc file formats: http://bioinformatics.ai.sri.com/ptools/flatfile-format.html
    * SBML FAQ: https://synonym.caltech.edu/documents/faq

    Args:
        neo4j: A :class:`.SBMLClient` instance.
        sbml: The path to the MetaCyc SBML file to convert.

    Attributes:
        db: A :class:`.SBMLClient` instance. This is connected to neo4j and used
        to perform all database operations. Should be closed after use.
        sbml_file: Filepath to the input SBML file.
    """

    # MIRIAM qualifiers (https://co.mbine.org/standards/qualifiers)
    # libsbml.BiolQualifierType_toString
    _bio_qualifiers = {
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

    def __init__(
        self,
        neo4j: SBMLClient,
        sbml: Union[str, Path],
    ):
        self.db = neo4j  # Neo4j driver
        self.sbml_file = validate_path(sbml)
        if not self.sbml_file:
            raise ValueError("Missing SBML file path.")

    def sbml_to_graph(self):
        """
        Populate Neo4j database with SBML data. The process is as follows:

        #. Parse the SBML file. All parsing errors are logged as warnings.
        #. Create the database and constraints.
        #. Feed the SBML file into the database. This will populate
           ``Compartment``, ``Reaction``, ``Compound``, ``GeneProduct``,
           ``GeneProductSet``, ``GeneProductComplex``, and ``RDF`` nodes.

        Nodes are created for each SBML element using ``MERGE`` statements:
        https://neo4j.com/docs/cypher-manual/current/clauses/merge/#merge-merge-with-on-create
        """
        # Read SBML file
        doc = self._read_sbml(self.sbml_file)
        model: libsbml.Model = doc.getModel()

        # Compartments
        compartments = self._collect_compartments(model.getListOfCompartments())
        self.db.create_nodes(
            "Compartment", compartments, self.db.default_cyphers["Compartment"]
        )

        # Compounds, i.e. metabolites, species
        compounds = self._collect_compounds(model.getListOfSpecies())
        self.db.create_nodes("Compound", compounds, self.db.default_cyphers["Compound"])

        # Gene products
        fbc: libsbml.FbcModelPlugin = model.getPlugin("fbc")
        gene_prods = self._collect_gene_products(fbc.getListOfGeneProducts())
        self.db.create_nodes(
            "GeneProduct", gene_prods, self.db.default_cyphers["GeneProduct"]
        )

        # Reactions
        rxns = model.getListOfReactions()
        reactions = self._collect_reactions(rxns)
        self.db.create_nodes(
            "Reaction",
            reactions,
            self.db.default_cyphers["Reaction"],
            progress_bar=True,
        )
        for r in tqdm(rxns, desc="Reaction-GeneProduct"):
            # Add associated gene products
            # This could be complicated where the child nodes could be:
            # 1. GeneProductRef
            # 2. fbc:or -> GeneProductRef
            # 3. fbc:and -> GeneProductRef
            # 4. Arbitrarily nested fbc:or / fbc:and within cases 2 and 3
            mcid = r.getMetaId()
            gpa: libsbml.GeneProductAssociation = r.getPlugin(
                "fbc"
            ).getGeneProductAssociation()
            if gpa is not None:
                node = gpa.getAssociation()
                self._add_sbml_gene_product_association_node(node, mcid)

    @staticmethod
    def _read_sbml(sbml_file: Path) -> libsbml.SBMLDocument:
        reader = libsbml.SBMLReader()
        doc = reader.readSBMLFromFile(sbml_file)
        logger.info("Finished reading SBML file")

        for i in range(doc.getNumErrors()):
            err = doc.getError(i)
            logger.warning(
                "SBML reader raised %s during parsing at line %d, column %d: %s",
                err.getSeverityAsString(),
                err.getLine(),
                err.getColumn(),
                err.getMessage().replace("\n", " "),
            )

        return doc

    @staticmethod
    def _collect_compartments(
        compartments: List[libsbml.Compartment],
    ) -> List[Dict[str, Union[str, Dict[str, str]]]]:
        nodes = [
            {"metaId": c.getMetaId(), "props": {"name": c.getName()}}
            for c in compartments
        ]

        return nodes

    def _collect_compounds(
        self,
        compounds: List[libsbml.Species],
    ) -> List[Dict[str, Union[str, Dict[str, str]]]]:
        nodes = []
        for c in tqdm(compounds, desc="Compounds"):
            node = {
                "metaId": c.getMetaId(),
                "props": {
                    "name": c.getName(),
                    "chemicalFormula": c.getPlugin("fbc").getChemicalFormula(),
                    "boundaryCondition": c.getBoundaryCondition(),
                    "hasOnlySubstanceUnits": c.getHasOnlySubstanceUnits(),
                    "constant": c.getConstant(),
                },
                "compartment": c.getCompartment(),
            }
            cvterms: List[libsbml.CVTerm] = c.getCVTerms()
            node["rdf"] = [self._cvterm_to_rdf(cvterm) for cvterm in cvterms]

            nodes.append(node)

        return nodes

    def _collect_gene_products(
        self, gene_prods: List[libsbml.GeneProduct]
    ) -> List[Dict[str, Union[str, Dict[str, str]]]]:
        nodes = []
        for gp in tqdm(gene_prods, desc="GeneProducts"):
            node = {
                "metaId": gp.getMetaId(),
                "props": {
                    "name": gp.getName(),
                    "label": gp.getLabel(),
                },
            }
            cvterms: List[libsbml.CVTerm] = gp.getCVTerms()
            node["rdf"] = [self._cvterm_to_rdf(cvterm) for cvterm in cvterms]

            nodes.append(node)

        return nodes

    def _collect_reactions(self, reactions: List[libsbml.Reaction]):
        nodes = []
        for r in tqdm(reactions, desc="Reactions"):
            r: libsbml.Reaction
            node = {
                "metaId": r.getMetaId(),
                "props": {
                    "name": r.getName(),
                    "fast": r.getFast(),
                },
            }
            cvterms: List[libsbml.CVTerm] = r.getCVTerms()
            node["rdf"] = [self._cvterm_to_rdf(cvterm) for cvterm in cvterms]

            # Add reactants and products
            reactants = r.getListOfReactants()
            node["reactants"] = self._link_reaction_to_compound(reactants)
            products = r.getListOfProducts()
            node["products"] = self._link_reaction_to_compound(products)

            nodes.append(node)

        return nodes

    def _cvterm_to_rdf(
        self, cvterm: libsbml.CVTerm
    ) -> Dict[str, Union[str, Dict[str, Union[str, List[str]]]]]:
        """Convert a CVTerm to RDF node.

        RDF in Annotation are in the form of triples:
        the model component to annotate (subject), the relationship between
        the model component and the annotation (predicate), and a term
        describing the component (object).

        Args:
            cvterm: The CVTerm that contains information to add to the RDF node.
        """
        # Get the biological qualifier type of the terms
        bio_qual = self._bio_qualifiers[cvterm.getBiologicalQualifierType()]
        # Get the content of each RDF term
        uris = [
            self._split_uri(cvterm.getResourceURI(i))
            for i in range(cvterm.getNumResources())
        ]

        props: Dict[str, Union[str, List[str]]] = {}
        for resource, identifier in uris:
            is_list_resource = resource == "ec-code"
            add_kv_to_dict(props, resource, identifier, as_list=is_list_resource)

        return {"bioQual": bio_qual, "rdf": props}

    def _link_reaction_to_compound(
        self,
        compounds: List[libsbml.SpeciesReference],
    ):
        """Link reactants or products to a reaction.

        Args:
            compounds: The list of compounds to link to the reaction.
        """
        cpds = [
            {
                "cpdId": cpd.getSpecies(),
                "props": {
                    "stoichiometry": cpd.getStoichiometry(),
                    "constant": cpd.getConstant(),
                },
            }
            for cpd in compounds
        ]
        return cpds

    def _add_sbml_gene_product_association_node(
        self,
        node: Union[libsbml.GeneProductRef, libsbml.FbcAnd, libsbml.FbcOr],
        source_id: str,
        source_label: str = "Reaction",
        edge_type: str = "hasGeneProduct",
        node_index: int = 0,
    ):
        """
        Add gene products to a reaction. When the added node is FbcAnd or FbcOr,
        recursively add the children. This means a custom ``metaId`` is
        constructed for the ``GeneProductComplex`` and ``GeneProductSet`` nodes
        corresponding to the ``FbcAnd`` and ``FbcOr`` nodes, respectively.

        Args:
            node: The GeneProductAssociation child node to add.
            source_id: The MetaCyc ID of the source node. This should be the
                MetaCyc ID of the ``Reaction`` node
            source_label: The label of the source node.
            edge_type: The type of edge to add. Should be one of
                ``hasGeneProduct``, ``hasComponent``, or ``hasMember``.
            node_index: The index of the current node. This is used to construct
                the ``metaId`` of the ``GeneProductComplex`` and
                ``GeneProductSet`` nodes.
        """
        # If there's no nested association, add the node directly
        if isinstance(node, libsbml.GeneProductRef):
            self.db.link_node_to_gene(
                source_label, source_id, node.getGeneProduct(), "GeneProduct", edge_type
            )

        # For nested associations, first add a `GeneProductComplex` or
        # `GeneProductSet` node, then recursively add the children
        elif isinstance(node, libsbml.FbcAnd):
            complex_id = f"{source_id}_complex{node_index}"
            self.db.link_node_to_gene(
                source_label, source_id, complex_id, "GeneProductComplex", edge_type
            )

            for i in range(node.getNumAssociations()):
                self._add_sbml_gene_product_association_node(
                    node.getAssociation(i),
                    complex_id,
                    "GeneProductComplex",
                    "hasComponent",
                    i,
                )
        elif isinstance(node, libsbml.FbcOr):
            eset_id = f"{source_id}_entityset{node_index}"
            self.db.link_node_to_gene(
                source_label, source_id, eset_id, "GeneProductSet", edge_type
            )

            for i in range(node.getNumAssociations()):
                self._add_sbml_gene_product_association_node(
                    node.getAssociation(i),
                    eset_id,
                    "GeneProductSet",
                    "hasMember",
                    i,
                )
        else:
            logging.error(f"Unhandled GeneProductAssociation type {type(node)}")
            raise ValueError

    @staticmethod
    def _split_uri(uri: str) -> Tuple[str, str]:
        """Split a URI into a namespace and an annotation term.

        Args:
            uri: URI to split.

        Returns:
            Tuple of namespace and annotation term.
        """
        # First three elements are from http://identifiers.org/
        res = uri.split("/")[3:]
        resource, identifier = res[0], res[1:]
        resource = resource.replace(".", "-")  # Ec-code
        identifier = "".join(identifier)

        # In some cases the identifier in the RDF nodes of the SBML file has a
        # prefix, e.g. META: or HUMAN: for BioCyc IDs. We don't want these as
        # they make later queries difficult.
        identifier = re.sub(r"^[A-Z]+:", "", identifier)
        return resource, identifier
