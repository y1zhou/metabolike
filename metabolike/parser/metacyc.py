import logging
import re
from itertools import groupby
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple, Union

import libsbml
import pandas as pd
from metabolike.db.metacyc import MetaDB
from tqdm import tqdm

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


class Metacyc:
    """
    Converting MetaCyc files to a Neo4j database.
    Documentation on the MetaCyc files and format FAQs can be found at:

    * MetaCyc data files download: https://metacyc.org/downloads.shtml
    * MetaCyc file formats: http://bioinformatics.ai.sri.com/ptools/flatfile-format.html
    * SBML FAQ: https://synonym.caltech.edu/documents/faq

    Args:
        neo4j: A :class:`.MetaDB` instance.
        sbml: The path to the MetaCyc SBML file to convert.
        reactions: The path to the ``reaction.dat`` file. If given,
            the file will be parsed and extra annotation on ``Reaction``
            nodes will be added.
        atom_mapping: The path to the ``atom-mappings-smiles.dat`` file. If
            given, the file will be parsed and chemical reactions in the ``SMILES``
            format will be added to the ``Reaction`` nodes.
        pathways: The path to the ``pathway.dat`` file. If given,
            the file will be parsed and pathway links will be added to the
            ``Reaction`` nodes.
        compounds: The path to the ``compound.dat`` file. If given, the
            file will be parsed and annotations on ``Compound`` nodes will
            be added.
        publications: The path to the ``publication.dat`` file. If given,
            the file will be parsed and annotations on ``Citation`` nodes
            will be added.
        classes: The path to the ``class.dat`` file. If given, the file
            will be parsed and annotations on ``Compartment``, ``Taxa``,
            and ``Compound`` nodes will be added.
        db_name: The name of the database to create. If None, the ``metaid``
            attribute of the SBML file is used.

    Attributes:
        db: A :class:`.MetaDB` instance. In :meth:`.setup`, this is connected
            to neo4j and used to perform all database operations. Should be
            closed after use.
        input_files: A dictionary of the paths to the input ``.dat`` files.
        missing_ids: A dictionary of sets of IDs that were not found in the
            input files. This is helpful for collecting IDs that appear to be
            in one class but are actually in another.
            :meth:`._report_missing_ids` can be used to print them out.
    """

    def __init__(
        self,
        neo4j: MetaDB,
        sbml: Union[str, Path],
        reactions: Optional[Union[str, Path]] = None,
        atom_mapping: Optional[Union[str, Path]] = None,
        pathways: Optional[Union[str, Path]] = None,
        compounds: Optional[Union[str, Path]] = None,
        publications: Optional[Union[str, Path]] = None,
        classes: Optional[Union[str, Path]] = None,
        db_name: Optional[str] = None,
    ):
        # Neo4j driver
        self.db = neo4j
        if db_name:
            self.db.use_database(db_name)

        # File paths
        self.input_files = {
            "sbml": validate_path(sbml),
            "reactions": validate_path(reactions),
            "atom_mapping": validate_path(atom_mapping),
            "pathways": validate_path(pathways),
            "compounds": validate_path(compounds),
            "publications": validate_path(publications),
            "classes": validate_path(classes),
        }
        logger.info(f"Input files: {self.input_files}")

        # Placeholder for missing IDs in the dat files
        self.missing_ids: Dict[str, Set[str]] = {
            "reactions": set(),
            "atom_mappings": set(),
            "pathways": set(),
            "compounds": set(),
            "publications": set(),
            "compartments": set(),
            "taxon": set(),
        }

    def setup(self, create_db: bool = True, **kwargs):
        """Enterpoint for setting up the database.

        The process is as follows:

        #. Parse the SBML file. All parsing errors are logged as warnings.
        #. If ``db_name`` is not given, use the ``metaid`` attribute of the SBML
           file to name the database.
        #. Create the database and constraints.
        #. Feed the SBML file into the database. This will populate
           ``Compartment``, ``Reaction``, ``Compound``, ``GeneProduct``,
           ``EntitySet``, ``Complex``, and ``RDF`` nodes.
        #. If ``reactions.dat`` is given, parse the file and add standard Gibbs
           free energy, standard reduction potential, reaction direction,
           reaction balance status, systematic name, comment attributes to
           ``Reaction`` nodes. Also link ``Reaction`` nodes to ``Pathway`` and
           ``Citation`` nodes.
           #. If ``pathways.dat`` is given:

              * Add synonyms, types, comments, common names to ``Pathway`` nodes.
              * Link ``Pathway`` nodes to super-pathway ``Pathway`` and taxonomy
                  ``Taxa`` nodes.
              * Link ``Reaction`` nodes within the pathway with
                  ``isPrecedingEvent`` relationships.
              * Link ``Pathway`` nodes with their rate limiting steps
                  (``Reaction`` nodes).
              * Link ``Pathway`` nodes with primary reactant and product
                  ``Compound`` nodes.
              * Link ``Pathway`` nodes with other ``Pathway`` nodes and annotate
                the shared ``Compound`` nodes.
              * For ``Reaction`` nodes, add ``isPrimary[Reactant|Product]InPathway``
                  labels in their links to ``Compound`` nodes.
        #. If ``atom-mappings-smiles.dat`` is given, parse the file and add
           SMILES_ mappings to ``Reaction`` nodes.

        #. If ``compounds.dat`` is given, parse the file and add standard Gibbs
           free energy, logP, molecular weight, monoisotopic molecular weight,
           polar surface area, pKa {1,2,3}, comment, and synonyms to ``Compound``
           nodes. Also add SMILES and INCHI strings to related ``RDF`` nodes, and
           link the ``Compound`` nodes to ``Citation`` nodes.
        #. If ``pubs.dat`` is given, parse the file and add DOI, PUBMED, MEDLINE IDs,
           title, source, year, URL, and ``REFERENT-FRAME`` to ``Citation`` nodes.
        #. If ``classes.dat`` is given, parse the file and:

           * Add common name and synonyms to ``Compartment`` nodes.
           * Add common name, strain name, comment, and synonyms to ``Taxa`` nodes.

        .. _SMILES: https://en.wikipedia.org/wiki/SMILES
        """
        # Read SBML file and set default database name
        doc = self._read_sbml(self.input_files["sbml"])
        model: libsbml.Model = doc.getModel()
        if not self.db.db_name:
            self.db.db_name = model.getMetaId().lower()

        # Setup Neo4j database and populate it with data
        self.db.use_database(self.db.db_name)
        self.db.setup_graph_db(create_db=create_db, **kwargs)
        self.sbml_to_graph(model)

        # Add additional information of reactions to the graph if given
        if self.input_files["reactions"] or self.input_files["atom_mapping"]:
            all_rxns = self.db.get_all_nodes("Reaction", "displayName")

            if self.input_files["reactions"]:
                logger.info("Adding additional reaction information to the graph")
                rxn_dat = self._read_dat_file(self.input_files["reactions"])

                for rxn in tqdm(all_rxns, desc="reactions.dat file"):
                    self.reaction_to_graph(rxn, rxn_dat)
                    logger.debug(f"Added extra info for reaction {rxn}")

                # Read pathways file if given
                if self.input_files["pathways"]:
                    logger.info("Creating pathway links")
                    pw_dat = self._read_dat_file(self.input_files["pathways"])

                    all_pws = self.db.get_all_nodes("Pathway", "mcId")
                    self.super_pathways = set()
                    for pw in tqdm(all_pws, desc="pathways.dat file"):
                        self.pathway_to_graph(pw, pw_dat, rxn_dat)
                        logger.debug(f"Added pathway annotation for {pw}")

                    # Annotate the super-pathways
                    self.super_pathways = self.super_pathways - set(all_pws)
                    while self.super_pathways:
                        spw = self.super_pathways.pop()
                        self.pathway_to_graph(spw, pw_dat, rxn_dat)
                        logger.debug(f"Added pathway annotation for {spw}")

                    del self.super_pathways

                    # Correct composite reactions
                    self.db.set_composite_reaction_labels()

            # SMILES reactions
            if self.input_files["atom_mapping"]:
                logger.info("Adding SMILES to reactions")
                smiles = pd.read_table(
                    self.input_files["atom_mapping"],
                    sep="\t",
                    header=None,
                    names=["rxn", "smiles"],
                )
                smiles: Dict[str, str] = smiles.set_index("rxn").to_dict()["smiles"]
                for rxn in tqdm(all_rxns, desc="atom_mapping.dat file"):
                    self.atom_mapping_to_graph(rxn, smiles)

        # Compounds in compounds.dat
        if self.input_files["compounds"]:
            logger.info("Annotating Compound nodes")
            cpd_dat = self._read_dat_file(self.input_files["compounds"])

            all_cpds = self.db.get_all_compounds()
            for cpd, biocyc in tqdm(all_cpds, desc="compounds.dat file"):
                self.compounds_to_graph(cpd, biocyc, cpd_dat)
                logger.debug(f"Added annotation for compound {cpd} with {biocyc}")

        # Read publications file if given
        if self.input_files["publications"]:
            logger.info("Annotating publications")
            pub_dat = self._read_dat_file(self.input_files["publications"])

            all_cits = self.db.get_all_nodes("Citation", "mcId")
            for cit in tqdm(all_cits, desc="pubs.dat file"):
                self.citation_to_graph(cit, pub_dat)
                logger.debug(f"Added annotation for citation {cit}")

        # Compartments, Taxon, and comments of Compounds in classes.dat
        if self.input_files["classes"]:
            logger.info("Adding common names from classes.dat")
            class_dat = self._read_dat_file(self.input_files["classes"])
            self.classes_to_graph(class_dat)

        self._report_missing_ids()

    def sbml_to_graph(self, model: libsbml.Model):
        """
        Populate Neo4j database with SBML data.

        Nodes are created for each SBML element using ``MERGE`` statements:
        https://neo4j.com/docs/cypher-manual/current/clauses/merge/#merge-merge-with-on-create
        """
        # Compartments
        logger.info("Creating Compartment nodes")
        for c in model.getListOfCompartments():
            c: libsbml.Compartment
            props = {"displayName": c.getName()}
            self.db.create_node("Compartment", c.getMetaId(), props)

        # Compounds, i.e. metabolites, species
        logger.info("Creating Compound nodes")
        for s in tqdm(model.getListOfSpecies(), desc="Compounds"):
            s: libsbml.Species
            # Basic properties
            mcid: str = s.getMetaId()
            props = {
                "displayName": s.getName(),
                "charge": s.getCharge(),
                "chemicalFormula": s.getPlugin("fbc").getChemicalFormula(),
                "boundaryCondition": s.getBoundaryCondition(),  # unchanged by reactions
                "hasOnlySubstanceUnits": s.getHasOnlySubstanceUnits(),
                "constant": s.getConstant(),
            }
            self.db.create_compound_node(s.getCompartment(), mcid, props)

            # Add RDF annotations
            logger.debug(f"Adding RDF nodes for Compound {mcid}")
            cvterms: List[libsbml.CVTerm] = s.getCVTerms()
            for cvterm in cvterms:
                self._add_sbml_rdf_node("Compound", mcid, cvterm)

        # Gene products
        logger.info("Creating GeneProduct nodes")
        fbc: libsbml.FbcModelPlugin = model.getPlugin("fbc")
        for gp in tqdm(fbc.getListOfGeneProducts(), desc="GeneProducts"):
            gp: libsbml.GeneProduct
            mcid = gp.getMetaId()
            props = {
                "displayName": gp.getName(),
                "label": gp.getLabel(),
            }
            self.db.create_node("GeneProduct", mcid, props)

            logger.debug(f"Adding RDF nodes for GeneProduct {mcid}")
            cvterms: List[libsbml.CVTerm] = gp.getCVTerms()
            for cvterm in cvterms:
                self._add_sbml_rdf_node("GeneProduct", mcid, cvterm)

        # Reactions
        logger.info("Creating Reaction nodes")
        for r in tqdm(model.getListOfReactions(), desc="Reactions"):
            r: libsbml.Reaction
            mcid = r.getMetaId()

            # Basic properties
            props = {
                "displayName": r.getName(),
                "reversible": r.getReversible(),
                "fast": r.getFast(),
            }
            self.db.create_node("Reaction", mcid, props)

            # Add RDF nodes
            logger.debug(f"Adding RDF nodes for Reaction {mcid}")
            cvterms: List[libsbml.CVTerm] = r.getCVTerms()
            for cvterm in cvterms:
                self._add_sbml_rdf_node("Reaction", mcid, cvterm)

            # Add reactants and products
            logger.debug(f"Adding reactants for Reaction {mcid}")
            reactants = r.getListOfReactants()
            self._link_reaction_to_compound(mcid, reactants, "Left")

            logger.debug(f"Adding products for Reaction {mcid}")
            products = r.getListOfProducts()
            self._link_reaction_to_compound(mcid, products, "Right")

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
            self._add_sbml_gene_product_association_node(node, mcid)

    def reaction_to_graph(self, rxn_id: str, rxn_dat: Dict[str, List[List[str]]]):
        """
        Parse one entry from the reaction attribute-value file, and add relevant
        information to the graph database in one transaction.

        Args:
            rxn_id: The *full* reaction ID from the graph database.
            rxn_dat: The reactions.dat file as an attribute-value list.
        """
        canonical_id = self._find_rxn_canonical_id(rxn_id, rxn_dat.keys())
        if canonical_id not in rxn_dat:
            self.missing_ids["reactions"].add(canonical_id)
            return
        lines = rxn_dat[canonical_id]
        props: Dict[str, Union[str, List[str]]] = {"canonicalId": canonical_id}

        for k, v in lines:
            # SYNONYMS is a special case because it is a list
            if k in {"SYNONYMS", "TYPES"}:
                _add_kv_to_dict(props, k, v, as_list=True)
            elif k in {
                "GIBBS-0",
                "STD-REDUCTION-POTENTIAL",
                "REACTION-DIRECTION",
                "REACTION-BALANCE-STATUS",
                "SYSTEMATIC-NAME",
                "COMMENT",  # TODO: link to other nodes
            }:
                _add_kv_to_dict(props, k, v, as_list=False)
            elif k == "IN-PATHWAY":
                self.db.link_reaction_to_pathway(rxn_id, v)
            elif k == "RXN-LOCATIONS":
                self.db.link_reaction_to_compartment(rxn_id, v)
            elif k == "CITATIONS":
                self.db.link_node_to_citation("Reaction", rxn_id, v)

        # Clean up props before writing to graph
        props = self._clean_props(
            props,
            num_fields=[
                _snake_to_camel(x) for x in ["GIBBS-0", "STD-REDUCTION-POTENTIAL"]
            ],
            enum_fields=[
                _snake_to_camel(x)
                for x in ["REACTION-BALANCE-STATUS", "REACTION-DIRECTION"]
            ],
        )
        self.db.add_props_to_node("Reaction", "displayName", rxn_id, props)

    def atom_mapping_to_graph(self, rxn_id: str, smiles: Dict[str, str]):
        """
        Args:
            rxn_id: The *full* reaction ID from the graph database.
            smiles: The reaction ID -> SMILES dictionary from the atom mapping file.
        """
        canonical_id = self._find_rxn_canonical_id(rxn_id, smiles.keys())
        if canonical_id not in smiles:
            self.missing_ids["atom_mappings"].add(canonical_id)
            return
        props = {"smiles_atom_mapping": smiles[canonical_id]}
        self.db.add_props_to_node("Reaction", "displayName", rxn_id, props)

    def pathway_to_graph(
        self,
        pw_id: str,
        pw_dat: Dict[str, List[List[str]]],
        rxn_dat: Dict[str, List[List[str]]],
    ):
        """
        Parse one entry from the pathway attribute-value file, and add relevant
        information to the graph database in one transaction.

        Args:
            pw_id: The pathway ID from the graph database.
            pw_dat: The pathways.dat file as an attribute-value list.
            rxn_dat: The reactions.dat file as an attribute-value list. This is
                used to find the annotation data for composite reactions marked
                as pathways.
        """
        if (pw_id not in pw_dat) and (pw_id not in rxn_dat):
            self.missing_ids["pathways"].add(pw_id)
            return

        # Deal with pathway nodes with reaction IDs (composite reactions)
        if pw_id in rxn_dat:
            self.reaction_to_graph(pw_id, rxn_dat)
            # Merge the Reaction node with the Pathway node under the same ID
            self.db.merge_nodes("Pathway", "Reaction", "mcId", "displayName", pw_id)
            return

        lines = pw_dat[pw_id]
        props: Dict[str, Union[str, List[str]]] = {"displayName": pw_id}
        for k, v in lines:
            # Pathway node properties
            if k in {"SYNONYMS", "TYPES"}:
                _add_kv_to_dict(props, k, v, as_list=True)
            elif k in {"COMMENT", "COMMON-NAME"}:
                _add_kv_to_dict(props, k, v, as_list=False)

            # Relationship with other nodes
            elif k in {"IN-PATHWAY", "SUPER-PATHWAYS"}:
                self.super_pathways.add(v)
                self.db.link_pathway_to_superpathway(pw_id, v)
            elif k == "CITATIONS":
                self.db.link_node_to_citation("Pathway", pw_id, v)
            elif k == "SPECIES":
                self.db.link_pathway_to_taxa(pw_id, v, "hasRelatedSpecies")
            elif k == "TAXONOMIC-RANGE":
                self.db.link_pathway_to_taxa(pw_id, v, "hasExpectedTaxonRange")
            elif k == "PREDECESSORS":
                rxns = v[2:-2].split('" "')  # leading (" and trailing ")
                r1, r2 = rxns[0], rxns[1:]  # could have multiple predecessors
                self.db.link_reaction_to_next_step(r1, r2, pw_id)
            elif k == "RATE-LIMITING-STEP":
                self.db.link_pathway_to_rate_limiting_step(pw_id, v)
            elif k in {"PRIMARY-PRODUCTS", "PRIMARY-REACTANTS"}:
                relationship = "Product" if k == "PRIMARY-PRODUCTS" else "Reactant"
                self.db.link_pathway_to_primary_compound(pw_id, v, relationship)
            elif k == "REACTION-LAYOUT":
                rxn_id, d = self._parse_reaction_layout(v)
                if not rxn_id:
                    continue
                # Two cypher statements for Reactant and Product compounds
                for side, compounds in d.items():
                    self.db.link_reaction_to_primary_compound(
                        rxn_id, compounds, pw_id, side
                    )
            elif k == "PATHWAY-LINKS":
                cpd_id, pw_ids, direction = self._parse_pathway_links(v)
                if pw_ids:
                    self.db.link_pathway_to_pathway(pw_id, pw_ids, direction, cpd_id)

        # Write Pathway node properties
        self.db.add_props_to_node("Pathway", "mcId", pw_id, props)

    def compounds_to_graph(
        self,
        cpd: str,
        biocyc: str,
        cpd_dat: Dict[str, List[List[str]]],
    ):
        """Annotate a compound node with data from the compound.dat file.

        Args:
            cpd: The ``displayName`` of the compound.
            biocyc: The biocyc id of the compound.
            cpd_dat: The compound.dat data.
            session: The neo4j session.
        """
        if biocyc not in cpd_dat:
            self.missing_ids["compounds"].add(biocyc)
            return
        lines = cpd_dat[biocyc]
        c_props: Dict[str, Union[str, List[str]]] = {}
        rdf_props: Dict[str, Union[str, List[str]]] = {}
        for k, v in lines:
            if k in {
                "COMMON-NAME",
                "GIBBS-0",
                "LOGP",
                "MOLECULAR-WEIGHT",
                "MONOISOTOPIC-MW",
                "POLAR-SURFACE-AREA",
                "PKA1",
                "PKA2",
                "PKA3",
                "COMMENT",
            }:
                _add_kv_to_dict(c_props, k, v, as_list=False)
            elif k == "SYNONYMS":
                _add_kv_to_dict(c_props, k, v, as_list=True)
            elif k in {"SMILES", "INCHI"}:
                _add_kv_to_dict(rdf_props, k, v, as_list=False)
            elif k == "DBLINKS":
                pass  # TODO: parse DBLINKS
            elif k == "CITATIONS":
                self.db.link_node_to_citation("Compound", cpd, v)

        c_props = self._clean_props(
            c_props,
            num_fields=[
                _snake_to_camel(x)
                for x in [
                    "GIBBS-0",
                    "LOGP",
                    "MOLECULAR-WEIGHT",
                    "MONOISOTOPIC-MW",
                    "POLAR-SURFACE-AREA",
                    "PKA1",
                    "PKA2",
                    "PKA3",
                ]
            ],
            enum_fields=[],
        )
        self.db.add_props_to_compound(cpd, c_props, rdf_props)

    def citation_to_graph(self, cit_id: str, pub_dat: Dict[str, List[List[str]]]):
        """Annotate a citation node with data from the publication.dat file.

        If there are multiple fields in the given ``cit_id``, then the fields are
        separated by colons. The first field is the citation ID, the second is the
        evidence type (in classes.dat), the third is not documented, and the
        fourth is the curator's name.

        In most cases the citation ID should match: ``PUB-[A-Z0-9]+$``,
        with a few exceptions containing double dashes, e.g. ``PUB--8``,
        or some dashes within author names, e.g. ``PUB-CHIH-CHING95``.

        Args:
            cit_id: The citation ``mcId`` property.
            pub_dat: The publication.dat data.
        """
        # TODO: deal with evidence frames. Evidence frames are in the form of
        # 10066805:EV-EXP-IDA:3354997111:hartmut
        pub_dat_id = re.sub(r"[\[\]\s,']", "", cit_id.split(":")[0].upper())
        pub_dat_id = re.sub(r"-(\d+)", r"\1", pub_dat_id)
        if pub_dat_id == "BOREJSZA-WYSOCKI94":
            pub_dat_id = "BOREJSZAWYSOCKI94"  # only exception with dash removed
        if not pub_dat_id:
            return
        pub_dat_id = "PUB-" + pub_dat_id
        if pub_dat_id not in pub_dat:
            self.missing_ids["publications"].add(pub_dat_id)
            return

        lines = pub_dat[pub_dat_id]
        props: Dict[str, Union[str, List[str]]] = {"citationId": pub_dat_id}
        for k, v in lines:
            # Pathway node properties
            if k == "AUTHORS":
                _add_kv_to_dict(props, k, v, as_list=True)
            elif k in {
                "DOI-ID",
                "PUBMED-ID",
                "MEDLINE-ID",
                "TITLE",
                "SOURCE",
                "YEAR",
                "URL",
                "REFERENT-FRAME",
            }:
                _add_kv_to_dict(props, k, v, as_list=False)

        self.db.add_props_to_node("Citation", "mcId", cit_id, props)

    def classes_to_graph(self, class_dat: Dict[str, List[List[str]]]):
        """Parse the classes.dat file.

        Add common name and synonyms to ``Compartment`` nodes. Add common name,
        strain name, comment, and synonyms to ``Taxa`` nodes.

        Args:
            class_dat: The classes.dat data.
        """
        # Common names for cell components
        all_cco = self.db.get_all_nodes("Compartment", "displayName")
        for cco in all_cco:
            if cco not in class_dat:
                self.missing_ids["compartments"].add(cco)
                continue
            props: Dict[str, Union[str, List[str]]] = {}
            for k, v in class_dat[cco]:
                if k == "COMMON-NAME":
                    _add_kv_to_dict(props, k, v, as_list=False)
                elif k == "SYNONYMS":
                    _add_kv_to_dict(props, k, v, as_list=True)

            self.db.add_props_to_node("Compartment", "displayName", cco, props)

        # Common names and synonyms for organisms. Some also have strain names
        all_taxon = self.db.get_all_nodes("Taxa", "mcId")
        for taxa in tqdm(all_taxon, desc="Taxon in classes.dat"):
            if taxa not in class_dat:
                self.missing_ids["taxon"].add(taxa)
                continue
            props: Dict[str, Union[str, List[str]]] = {}
            for k, v in class_dat[taxa]:
                if k in {"COMMON-NAME", "STRAIN-NAME", "COMMENT"}:
                    _add_kv_to_dict(props, k, v, as_list=False)
                elif k == "SYNONYMS":
                    _add_kv_to_dict(props, k, v, as_list=True)
                # TODO: TYPES links the taxon

            self.db.add_props_to_node("Taxa", "mcId", taxa, props)

        # TODO: Evidence code in citations are in the `Evidence` attr

    @staticmethod
    def _read_sbml(sbml_file: Path) -> libsbml.SBMLDocument:
        reader = libsbml.SBMLReader()
        metacyc = reader.readSBMLFromFile(sbml_file)
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

    def _add_sbml_rdf_node(self, node_label: str, mcid: str, cvterm: libsbml.CVTerm):
        """Create RDF node and link it to the given SBML node.

        RDF in Annotation are in the form of triples:
        the model component to annotate (subject), the relationship between
        the model component and the annotation (predicate), and a term
        describing the component (object).

        Args:
            node_label: The label of the SBML node to annotate. Should be one
                of the labels defined in ``NODE_LABELS``.
            mcid: The MetaCyc ID of the SBML node to annotate.
            cvterm: The CVTerm that contains information to add to the RDF node.
        """
        # Get the biological qualifier type of the terms
        bio_qual = BIO_QUALIFIERS[cvterm.getBiologicalQualifierType()]
        # Get the content of each RDF term
        uris = [
            self._split_uri(cvterm.getResourceURI(i))
            for i in range(cvterm.getNumResources())
        ]

        props: Dict[str, Union[str, List[str]]] = {}
        for resource, identifier in uris:
            is_list_resource = resource == "ec-code"
            _add_kv_to_dict(props, resource, identifier, as_list=is_list_resource)

        self.db.link_node_to_rdf(node_label, mcid, bio_qual, props)

    def _link_reaction_to_compound(
        self,
        reaction_id: str,
        compounds: List[libsbml.SpeciesReference],
        compound_type: str,
    ):
        """Link reactants or products to a reaction.

        Args:
            reaction_id: The MetaCyc ID of the reaction.
            compounds: The list of compounds to link to the reaction.
            compound_type: The type of compound to link to the reaction. Should
                be one of "Left" or "Right".
        """
        if compound_type not in ["Left", "Right"]:
            raise ValueError(f"Invalid compound type: {compound_type}")
        for cpd in compounds:
            props = {
                "stoichiometry": cpd.getStoichiometry(),
                "constant": cpd.getConstant(),
            }
            self.db.link_reaction_to_compound(
                reaction_id, cpd.getSpecies(), compound_type, props
            )

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
        recursively add the children. This means a custom ``mcId`` is constructed
        for the ``Complex`` and ``EntitySet`` nodes corresponding to the ``FbcAnd``
        and ``FbcOr`` nodes, respectively.

        Args:
            node: The GeneProductAssociation child node to add.
            source_id: The MetaCyc ID of the source node. This should be the
                MetaCyc ID of the ``Reaction`` node
            source_label: The label of the source node.
            edge_type: The type of edge to add. Should be one of ``hasGeneProduct``,
                ``hasComponent``, or ``hasMember``.
            node_index: The index of the current node. This is used to construct
                the ``mcId`` of the ``Complex`` and ``EntitySet`` nodes.
        """
        # If there's no nested association, add the node directly
        if isinstance(node, libsbml.GeneProductRef):
            self.db.link_node_to_gene(
                source_label, source_id, node.getGeneProduct(), "GeneProduct", edge_type
            )

        # For nested associations, first add a `Complex` or `EntitySet` node,
        # then recursively add the children
        elif isinstance(node, libsbml.FbcAnd):
            complex_id = f"{source_id}_complex{node_index}"
            self.db.link_node_to_gene(
                source_label, source_id, complex_id, "Complex", edge_type
            )

            for i in range(node.getNumAssociations()):
                self._add_sbml_gene_product_association_node(
                    node.getAssociation(i),
                    complex_id,
                    "Complex",
                    "hasComponent",
                    i,
                )
        elif isinstance(node, libsbml.FbcOr):
            eset_id = f"{source_id}_entityset{node_index}"
            self.db.link_node_to_gene(
                source_label, source_id, eset_id, "EntitySet", edge_type
            )

            for i in range(node.getNumAssociations()):
                self._add_sbml_gene_product_association_node(
                    node.getAssociation(i),
                    eset_id,
                    "EntitySet",
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

    @staticmethod
    def _read_dat_file(filepath: Union[str, Path]) -> Dict[str, List[List[str]]]:
        # `rb` to bypass the 0xa9 character
        with open(filepath, "rb") as f:
            lines = [l.decode("utf-8", "ignore").strip() for l in f]
            # Also remove empty lines and comments on the top of the file
            lines = [l for l in lines if l and (not l.startswith("#"))]

        # Split entries based on `//`
        docs: Dict[str, List[List[str]]] = {}
        for k, g in groupby(lines, key=lambda x: x != "//"):
            if not k:
                continue

            doc = list(g)
            # Concatenate attributes with multiple lines (mostly COMMENT)
            doc_txt = "\n".join(doc)
            doc_txt = doc_txt.replace("\n/", " ")
            doc = doc_txt.split("\n")

            # Split key-attribute pairs
            doc = [l.split(" - ", maxsplit=1) for l in doc]
            uniq_id = doc[0][1]

            # Remove empty values like "SYNONYMS - " in ORG-6026
            doc = [x for x in doc[1:] if len(x) == 2]

            docs[uniq_id] = doc

        return docs

    @staticmethod
    def _find_rxn_canonical_id(rxn_id: str, all_ids: Iterable[str]) -> str:
        """Find the canonical ID for a reaction.

        Some reactions have longer ID forms in ``metabolic-reactions.xml`` than
        in ``reactions.dat`` or ``pathways.dat``. For example,
        ``F16BDEPHOS-RXN`` has two counterparts in ``metabolic-reactions.xml``:
        ``F16BDEPHOS-RXN[CCO-PERI-BAC]-FRUCTOSE-16-DIPHOSPHATE/WATER//FRUCTOSE-6P/Pi.60.``
        and
        ``F16BDEPHOS-RXN[CCO-CYTOSOL]-FRUCTOSE-16-DIPHOSPHATE/WATER//FRUCTOSE-6P/Pi.59.``

        This helper function extracts the leading part from the full ID.

        Args:
            rxn_id: The MetaCyc ID of the reaction.
            all_ids: All UNIQUE-IDs in the reactions.dat file.

        Returns:
            The canonical ID for the reaction.
        """
        # See if the rxn_id is in reaction.dat
        if rxn_id in all_ids:
            return rxn_id

        # If not, extract the ID from the leading part of rxn_id
        match_canonical_id = re.match(
            r"((TRANS-)?(RXN[A-Z\d]*)-\d+)|([A-Z\d.\-\+]+RXN)", rxn_id
        )
        if match_canonical_id:
            canonical_id = match_canonical_id.group(0)
            return canonical_id
        else:
            raise ValueError(f"rxn_id has no canonical form: {rxn_id}")

    @staticmethod
    def _clean_props(
        props: Dict[str, Any], num_fields: Iterable[str], enum_fields: Iterable[str]
    ) -> Dict[str, Any]:
        """Normalize properties to be used in Cypher.

        Args:
            props: Properties to normalize.
            num_fields: Fields that should be converted to float numbers.
            enum_fields: Fields that should be converted to alphanumerical strings.

        Returns:
            A dictionary with normalized properties.
        """
        for f in num_fields:
            if f in props:
                props[f] = float(props[f])

        enum_pattern = re.compile(r"\W+")
        for f in enum_fields:
            if f in props:
                props[f] = props[f].replace("-", "_").lower()
                props[f] = enum_pattern.sub("", props[f])

        return props

    @staticmethod
    def _parse_reaction_layout(s: str) -> Tuple[str, Dict[str, List[str]]]:
        """Parse the reaction layout from the ``reactions.dat`` file.

        The attribute has the following format: ``(<rxn-id> (:LEFT-PRIMARIES
        <cpd-id> ...) (:DIRECTION :[L2R|R2L]) (:RIGHT-PRIMARIES <cpd-id> ...))``

        Args:
            s: The layout string.

        Returns:
            The reaction canonical ID, and a dictionary with the parsed layout
            containing the following keys: - ``Reactant``: List of primary
            reactant compounds. - ``Product``: List of primary product
            compounds.
        """
        m = re.compile(
            r"\(([^\(]+) "  # reaction ID
            r"\(:LEFT-PRIMARIES ([^\)]+)\) "  # left primary compounds
            r"\(:DIRECTION :(L2R|R2L)\) "  # reaction direction
            r"\(:RIGHT-PRIMARIES ([^\)]+)\)\)"
        )
        res = m.fullmatch(s)
        if not res:
            # Sometimes there's no left/right primaries and only a direction
            # when the rxn_id is actually a pathway ID
            return "", {}

        rxn_id, l_prim, direction, r_prim = res.groups()
        # Use direction to determine which side are the reactants
        if direction == "L2R":
            reactants = l_prim.split(" ")
            products = r_prim.split(" ")
        else:
            reactants = r_prim.split(" ")
            products = l_prim.split(" ")

        d = {
            "Reactant": reactants,
            "Product": products,
        }

        return rxn_id, d

    @staticmethod
    def _parse_pathway_links(s: str) -> Tuple[str, List[str], str]:
        """
        Parse ``PATHWAY-LINKS`` from the ``pathways.dat`` file. This connects
        the ``Pathway`` nodes in the graph through shared ``Compound`` nodes.
        The format of the input string has several possibilities. The simplest
        is when the pathway points to one or more target pathways:

            (<compound-id> <target-pathway-ids ...>)

        These ``target-pathway-ids`` are inheriently *outgoing* links. Sometimes
        the direction is explicitly specified. When there are multiple pathway
        IDs and only one is specified as ``OUTGOING``, all the other pathways
        are also *outgoing* links.

            (<compound-id> (<target-pathway-id> . :OUTGOING))
            (<compound-id> (<target-pathway-id> . :INCOMING))

        The compound or pathway ID could also be a frame ID instead of the unique
        ID of the object, in which case the ID is wrapped in pipe symbols. This
        mostly happens for compound IDs, and when it happens for pathway IDs it's
        often a pathway type instead of the unique ID.

            (|<compound-frame-id>| <target-pathway-id>)

        Finally, the target pathway ID could be some random string wrapped in
        quotes. We can safely ignore most of these (e.g. "dietary input",
        "release from the lysosome"), but in several cases valid pathway IDs
        are also wrapped in quotes.

        Args:
            s: The pathway links string to parse.

        Returns:
            A tuple with the first element being the compound ID, the second
            element being a list of pathway IDs, and the third element being
            the direction of the link. The list of pathway IDs could be empty
            because non-canonical pathway IDs are dropped.
        """
        # Extract the compound ID at the beginning of the string
        cpd_id_rgx = re.compile(r"^\(([^\s]+) ")
        res = cpd_id_rgx.match(s)
        if not res:
            logging.warning(f"Pathway links string doesn't have a compound: {s}")
            return "", [], ""
        cpd = res.group(1)
        s = s[res.end() : -1]  # remove cpd ID and closing parenthesis

        # Extract pathway links recursively
        pathways: List[str] = []
        direction = ""
        pw_id = r"\"[^\"]+\"|[^\s]+"
        directed_pw_rgx = re.compile(rf"\(({pw_id}) \. :(INCOMING|OUTGOING)\)")
        pw_rgx = re.compile(pw_id)

        while s:
            # Directed pathway links wrapped in parentheses
            if s[0] == "(":
                m = directed_pw_rgx.match(s)
                pathways.append(m.group(1))
                direction = m.group(2)
            # Regular pathway links separated by spaces
            else:
                m = pw_rgx.match(s)
                pathways.append(m.group(0))
            s = s[m.end() :].strip()

        # Final cleanup
        cpd = cpd.replace("|", "")
        pathways = [
            p.replace('"', "") for p in pathways if (" " not in p) and p.isupper()
        ]  # remove frame IDs and quoted annotations
        direction = "INCOMING" if direction == "INCOMING" else "OUTGOING"
        return cpd, pathways, direction

    def _report_missing_ids(self):
        for datfile, ids in self.missing_ids.items():
            if ids:
                logger.warning(f"The following IDs were not found in {datfile}: {ids}")


def validate_path(filepath: Optional[Union[str, Path]]) -> Optional[Path]:
    if not filepath:
        return None
    f = Path(filepath).expanduser().resolve()
    if not f.is_file():
        logger.error(f"File does not exist: {f}")
        raise FileNotFoundError(str(f))

    return f


def _snake_to_camel(s: str, sep: str = "-") -> str:
    """Convert snake_case to CamelCase.

    Args:
        s: String to convert.
        sep: Separator to use.

    Returns:
        camelCase string.
    """
    s_camel = [x.capitalize() for x in s.split(sep)]
    s_camel[0] = s_camel[0].lower()
    return "".join(s_camel)


def _add_kv_to_dict(
    d: Dict[str, Any], k: str, v: Any, as_list: bool = False
) -> Dict[str, Any]:
    """Add a key-value pair to a dictionary.

    Args:
        d: Dictionary to add the key-value pair to.
        k: Key to add.
        v: Value to add.

    Returns:
        Dictionary with the key-value pair added. If the key already exists,
        the value is saved as a list.
    """
    k_camel = _snake_to_camel(k)
    if as_list:
        d.setdefault(k_camel, []).append(v)
    else:
        d[k_camel] = v

    return d
