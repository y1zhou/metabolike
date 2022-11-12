import logging
import re
from collections.abc import Iterable
from pathlib import Path
from typing import Union

import libsbml
from tqdm import tqdm

from metabolike.utils import add_kv_to_dict, validate_path

logger = logging.getLogger(__name__)


class SBMLParser:
    """Converting MetaCyc files to a Neo4j database. Documentation on the MetaCyc files and format
    FAQs can be found at:

    * [MetaCyc data files download](https://metacyc.org/downloads.shtml)
    * [MetaCyc file formats](https://bioinformatics.ai.sri.com/ptools/flatfile-format.html)
    * [SBML FAQ](https://synonym.caltech.edu/documents/faq)

    Args:
        sbml: The path to the MetaCyc SBML file to convert.

    Attributes:
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

    def __init__(self, sbml: Union[str, Path]):
        self.sbml_file = validate_path(sbml)
        if not self.sbml_file:
            raise ValueError("Missing SBML file path.")
        logger.info(f"SBML file: {self.sbml_file}")

    @staticmethod
    def read_sbml(sbml_file: Path) -> libsbml.SBMLDocument:
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
    def collect_compartments(
        compartments: Iterable[libsbml.Compartment],
    ) -> list[dict[str, Union[str, dict[str, str]]]]:
        return [{"metaId": c.getId(), "props": {"name": c.getName()}} for c in compartments]

    def collect_compounds(
        self,
        compounds: Iterable[libsbml.Species],
    ) -> list[dict[str, Union[str, dict[str, str]]]]:
        nodes = {}
        for c in tqdm(compounds, desc="Compounds"):
            if (metaId := c.getId()) in nodes:
                continue
            node = {
                "metaId": metaId,
                "props": {
                    "name": c.getName(),
                    "boundaryCondition": c.getBoundaryCondition(),
                    "hasOnlySubstanceUnits": c.getHasOnlySubstanceUnits(),
                    "constant": c.getConstant(),
                },
                "compartment": c.getCompartment(),
            }
            if (fbc := c.getPlugin("fbc")) is not None:
                node["props"]["chemicalFormula"] = fbc.getChemicalFormula()

            cvterms: list[libsbml.CVTerm] = c.getCVTerms()
            node["rdf"] = [self._cvterm_to_rdf(cvterm) for cvterm in cvterms]

            nodes[metaId] = node

        return list(nodes.values())

    def collect_gene_products(
        self, gene_prods: Iterable[libsbml.GeneProduct]
    ) -> list[dict[str, Union[str, dict[str, str]]]]:
        nodes = {}
        for gp in tqdm(gene_prods, desc="GeneProducts"):
            if (metaId := gp.getId()) in nodes:
                continue
            node = {
                "metaId": metaId,
                "props": {
                    "name": gp.getName(),
                    "label": gp.getLabel(),
                },
            }
            cvterms: list[libsbml.CVTerm] = gp.getCVTerms()
            node["rdf"] = [self._cvterm_to_rdf(cvterm) for cvterm in cvterms]

            nodes[metaId] = node

        return list(nodes.values())

    def collect_reactions(
        self, reactions: Iterable[libsbml.Reaction]
    ) -> list[dict[str, Union[str, dict[str, str]]]]:
        nodes = {}
        for r in tqdm(reactions, desc="Reactions"):
            r: libsbml.Reaction
            if (metaId := r.getId()) in nodes:
                continue
            node = {
                "metaId": metaId,
                "props": {
                    "name": r.getName(),
                    "fast": r.getFast(),
                    "reversible": r.getReversible(),
                },
            }
            cvterms: list[libsbml.CVTerm] = r.getCVTerms()
            node["rdf"] = [self._cvterm_to_rdf(cvterm) for cvterm in cvterms]

            # Add reactants and products
            reactants = r.getListOfReactants()
            node["reactants"] = self._get_reaction_compounds(reactants)
            products = r.getListOfProducts()
            node["products"] = self._get_reaction_compounds(products)

            nodes[metaId] = node

        return list(nodes.values())

    @staticmethod
    def collect_reverse_reactions(
        reactions: Iterable[libsbml.Reaction],
    ) -> list[dict[str, Union[str, dict[str, str]]]]:
        rev_reactions = []
        for r in reactions:
            if r["props"].get("reversible"):
                node = {
                    "metaId": f"rev-{r['metaId']}",
                    "props": {
                        "name": f"rev-{r['props']['name']}",
                        "fast": r["props"]["fast"],
                    },
                    "reactants": r["products"],
                    "products": r["reactants"],
                }
                rev_reactions.append(node)

        return rev_reactions

    @staticmethod
    def collect_groups(
        groups: Iterable[libsbml.Group],
    ) -> list[dict[str, Union[str, dict[str, str], list[str]]]]:
        return [
            {
                "metaId": g.getId(),
                "props": {
                    "name": g.getName(),
                    "kind": g.getKindAsString(),
                },
                "members": [m.getIdRef() for m in g.getListOfMembers()],
            }
            for g in groups
        ]

    def _cvterm_to_rdf(
        self, cvterm: libsbml.CVTerm
    ) -> dict[str, Union[str, dict[str, Union[str, list[str]]]]]:
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
        uris = [self._split_uri(cvterm.getResourceURI(i)) for i in range(cvterm.getNumResources())]

        props: dict[str, Union[str, list[str]]] = {}
        for resource, identifier in uris:
            is_list_resource = resource == "ec-code"
            add_kv_to_dict(props, resource, identifier, as_list=is_list_resource)

        return {"bioQual": bio_qual, "rdf": props}

    @staticmethod
    def _get_reaction_compounds(
        compounds: list[libsbml.SpeciesReference],
    ) -> list[dict[str, Union[str, dict[str, str]]]]:
        """Link reactants or products to a reaction.

        Args:
            compounds: The list of compounds to link to the reaction.
        """
        return [
            {
                "cpdId": cpd.getSpecies(),
                "props": {
                    "stoichiometry": cpd.getStoichiometry(),
                    "constant": cpd.getConstant(),
                },
            }
            for cpd in compounds
        ]

    def collect_reaction_gene_product_links(self, reactions: Iterable[libsbml.Reaction]):
        """Add gene products to a reaction. This could be complicated where the child nodes could
        be:

        #. GeneProductRef
        #. fbc:or -> GeneProductRef
        #. fbc:and -> GeneProductRef
        #. fbc:or -> {fbc:and -> GeneProductRef, GeneProductRef}

        Args:
            reactions: An iterable of SBML reactions.
        """
        # We need to know the metaId of children of each geneSet / geneComplex.
        # Also need to name the geneSet / geneComplex nodes.
        # Note that the same set/complex can catalyze multiple reactions.
        gene_sets: dict[str, set[str]] = {}
        gene_complexes: dict[str, set[str]] = {}

        # For each reaction, we store its associated gene products
        reaction_genes = {}
        for r in tqdm(reactions, desc="Associated GeneProduct"):
            mcid = r.getId()
            gpa = r.getPlugin("fbc").getGeneProductAssociation()
            if gpa is not None:
                gpa: libsbml.GeneProductAssociation
                node = gpa.getAssociation()
                reaction_genes[mcid] = self._parse_gene_association_nodes(
                    node, gene_sets, gene_complexes
                )

        return reaction_genes, gene_sets, gene_complexes

    def _parse_gene_association_nodes(
        self,
        node: Union[libsbml.GeneProductRef, libsbml.FbcAnd, libsbml.FbcOr],
        gene_sets: dict[str, set[str]],
        gene_complexes: dict[str, set[str]],
    ):
        """Get unique gene product sets/complexes from a GeneProductAssociation.

        Returns:
            The metaId of the associated gene product node of the reaction.
        """
        if isinstance(node, libsbml.GeneProductRef):
            # Add the gene product to the reaction
            reaction_link_id = node.getGeneProduct()
        elif isinstance(node, libsbml.FbcAnd):
            reaction_link_id = self._parse_gene_product_complex(node, gene_sets, gene_complexes)
        elif isinstance(node, libsbml.FbcOr):
            reaction_link_id = self._parse_gene_product_set(node, gene_sets, gene_complexes)
        else:
            raise ValueError(f"Unhandled GeneProductAssociation type {type(node)}")

        return reaction_link_id

    def _parse_gene_product_complex(
        self,
        node: libsbml.FbcAnd,
        gene_sets: dict[str, set[str]],
        gene_complexes: dict[str, set[str]],
    ):
        components = set()
        for i in range(node.getNumAssociations()):
            g = node.getAssociation(i)
            if isinstance(g, libsbml.GeneProductRef):
                components.add(g.getGeneProduct())
            elif isinstance(g, libsbml.FbcOr):
                gene_set = self._parse_gene_product_set(g, gene_sets, gene_complexes)
                components.add(gene_set)
            else:
                raise ValueError(node.getId())
        if components:
            return self._add_gene_product_group(components, gene_complexes, "GeneProductComplex")

    def _parse_gene_product_set(
        self,
        node: libsbml.FbcOr,
        gene_sets: dict[str, set[str]],
        gene_complexes: dict[str, set[str]],
    ):
        members = set()
        for i in range(node.getNumAssociations()):
            member = node.getAssociation(i)
            if isinstance(member, libsbml.GeneProductRef):
                members.add(member.getGeneProduct())
            elif isinstance(member, libsbml.FbcAnd):
                gene_complex = self._parse_gene_product_complex(member, gene_sets, gene_complexes)
                members.add(gene_complex)
            else:
                raise ValueError(node.getId())
        if members:
            return self._add_gene_product_group(members, gene_sets, "GeneProductSet")

    @staticmethod
    def _add_gene_product_group(ids: set[str], gene_group: dict[str, set[str]], group_type: str):
        """Add a set of genes to a dict of gene sets / complexes."""
        ids.discard(None)
        for k, v in gene_group.items():
            if ids == v:
                return k

        new_k = f"{group_type}{len(gene_group)}"
        gene_group[new_k] = ids
        return new_k

    @staticmethod
    def _split_uri(uri: str) -> tuple[str, str]:
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
