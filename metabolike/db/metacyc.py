import logging
from typing import Any, Dict, List, Optional, Tuple

from metabolike.parser import MetacycParser
from tqdm import tqdm

from .sbml import SBMLClient

logger = logging.getLogger(__name__)

_reactions_dat_cypher = """
UNWIND $batch_nodes AS n
  MERGE (r:Reaction {name: n.name})
    ON CREATE SET r += n.props
    ON MATCH SET r += n.props
  FOREACH (pw IN n.inPathway |
    MERGE (p:Pathway {metaId: pw})
      ON CREATE SET p.name = pw
    MERGE (p)-[:hasReaction]->(r)
  )
  FOREACH (cpt IN n.rxnLocations |
    MERGE (c:Compartment {metaId: cpt})
    MERGE (r)-[:hasCompartment]->(c)
  )
  FOREACH (cit IN n.citations |
    MERGE (c:Citation {metaId: cit})
    MERGE (r)-[:hasCitation]->(c)
  );
"""

_compounds_dat_cypher = """
UNWIND $batch_nodes AS n
  MATCH (c:Compound {name: n.name})
    SET c += n.props
  FOREACH (cit IN n.citations |
    MERGE (x:Citation {metaId: cit})
    MERGE (c)-[:hasCitation]->(x)
  );
"""

_pathways_dat_cypher = """
UNWIND $batch_nodes AS n
  MERGE (pw:Pathway {metaId: n.metaId})
    ON CREATE SET pw.name = n.metaId, pw += n.props
    ON MATCH SET pw += n.props
  FOREACH (cit IN n.citations |
    MERGE (c:Citation {metaId: cit})
    MERGE (pw)-[:hasCitation]->(c)
  )
  FOREACH (taxa IN n.species |
    MERGE (t:Taxa {metaId: taxa})
    MERGE (pw)-[:hasRelatedSpecies]->(t)
  )
  FOREACH (taxa IN n.taxonomicRange |
    MERGE (t:Taxa {metaId: taxa})
    MERGE (pw)-[:hasExpectedTaxonRange]->(t)
  )
  FOREACH (rxn IN n.rateLimitingStep |
    MERGE (pw)-[l:hasReaction]->(:Reaction {name: rxn})
      ON MATCH SET l.isRateLimitingStep = true
  )
  FOREACH (cpd IN n.primaryReactants |
    MERGE (c:Compound)-[:hasRDF {bioQualifier: 'is'}]->(:RDF {biocyc: cpd})
    MERGE (pw)-[:hasPrimaryReactant]->(c)
  )
  FOREACH (cpd IN n.primaryProducts |
    MERGE (c:Compound)-[:hasRDF {bioQualifier: 'is'}]->(:RDF {biocyc: cpd})
    MERGE (pw)-[:hasPrimaryProduct]->(c)
  )
  FOREACH (super_pw IN n.superPathways |
    MERGE (spw:Pathway {metaId: super_pw})
      ON CREATE SET spw.name = super_pw
    MERGE (spw)-[:hasSubPathway]->(pw)
  )
  FOREACH (super_pw IN n.inPathway |
    MERGE (spw:Pathway {metaId: super_pw})
      ON CREATE SET spw.name = super_pw
    MERGE (spw)-[:hasSubPathway]->(pw)
  )
  FOREACH (pred IN n.predecessors |
    MERGE (r1: Reaction {canonicalId: pred.r1})
    FOREACH (rxn IN pred.r2 |
      MERGE (r2:Reaction {canonicalId: rxn})
      MERGE (r2)-[l:isPrecedingEvent]->(r1)
        ON CREATE SET l.hasRelatedPathway = n.metaId
    )
  );
"""


class MetacycClient(SBMLClient):
    metacyc_default_cyphers = {
        "reactions": _reactions_dat_cypher,
        "compounds": _compounds_dat_cypher,
        "pathways": _pathways_dat_cypher,
    }

    def metacyc_to_graph(self, parser: MetacycParser):
        """Enterpoint for setting up the database.

        The process is as follows:

        #. Parse the SBML file. All parsing errors are logged as warnings.
        #. If ``db_name`` is not given, use the ``metaid`` attribute of the SBML
           file to name the database.
        #. Create the database and constraints.
        #. Feed the SBML file into the database. This will populate
           ``Compartment``, ``Reaction``, ``Compound``, ``GeneProduct``,
           ``GeneProductSet``, ``GeneProductComplex``, and ``RDF`` nodes.
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
        # Populate Neo4j database with data from SBML file
        self.sbml_to_graph(parser)

        rxn_dat = self._reactions_dat_to_graph(parser)
        if rxn_dat is not None:
            self._pathways_to_graph(rxn_dat, parser)

        self._smiles_dat_to_graph(parser)
        self._compounds_dat_to_graph(parser)
        self._citations_dat_to_graph(parser)
        self._classes_dat_to_graph(parser)

        parser.report_missing_ids()

    def _reactions_dat_to_graph(
        self, parser: MetacycParser
    ) -> Optional[Dict[str, List[List[str]]]]:
        """
        Parse the ``reactions.dat`` file and:

        * Add properties to ``Reaction`` nodes.
        * Link ``Reaction`` nodes to ``Pathway``, ``Citation``, and
        ``Compartment`` nodes.

        Returns:
            A dictionary of the parsed file. This is useful for passing into
            :meth:`pathways_to_graph` if the ``pathways.dat`` file is also
            given.
        """
        if not parser.input_files["reactions"]:
            logger.warning("No reactions.dat file given")
            return
        logger.info(f"Parsing {parser.input_files['reactions']}")
        all_rxns = self.get_all_nodes("Reaction", "name")
        rxn_dat = parser.read_dat_file(parser.input_files["reactions"])

        rxn_nodes = parser.collect_reactions_dat_nodes(all_rxns, rxn_dat)
        self._set_metaid_constraints("Pathway")
        self._set_metaid_constraints("Citation")
        self.create_nodes(
            "reactions.dat",
            rxn_nodes,
            self.metacyc_default_cyphers["reactions"],
            batch_size=100,
            progress_bar=True,
        )

        # Fix the direction of reactions
        self._fix_reaction_direction()
        return rxn_dat

    def _pathways_to_graph(
        self, rxn_dat: Dict[str, List[List[str]]], parser: MetacycParser
    ):
        """Parse the ``pathways.dat`` file.

        For details, see :meth:`setup`.
        """
        if not parser.input_files["pathways"]:
            logger.warning("No pathways.dat file given")
            return

        logger.info(f"Parsing {parser.input_files['pathways']}")
        pw_dat = parser.read_dat_file(parser.input_files["pathways"])
        all_pws = self.get_all_nodes("Pathway", "metaId")

        pw_nodes, comp_rxn_nodes = parser.collect_pathways_dat_nodes(
            all_pws, pw_dat, rxn_dat
        )
        for n in pw_nodes:
            if "species" in n or "taxonomicRange" in n:
                self._set_metaid_constraints("Taxa")
                break
        self._fix_composite_reaction_nodes(comp_rxn_nodes)

        # Annotate regular pathway nodes
        all_rxns = set(self.get_all_nodes("Reaction", "name"))
        pw_nodes = parser.fix_pathway_nodes(pw_nodes, all_rxns)
        self.create_nodes(
            "Pathway",
            pw_nodes,
            self.metacyc_default_cyphers["pathways"],
            batch_size=100,
            progress_bar=True,
        )

        for node in tqdm(pw_nodes, desc="Annotating pathways"):
            pw_id = node["metaId"]
            if rxn_layout := node.get("reactionLayout"):
                for v in rxn_layout:
                    # Two cypher statements for Reactant and Product compounds
                    self._link_reaction_to_primary_compound(
                        v["reaction"], v["reactants"], pw_id, "Reactant"
                    )
                    self._link_reaction_to_primary_compound(
                        v["reaction"], v["products"], pw_id, "Product"
                    )

            if pw_links := node.get("pathwayLinks"):
                for v in pw_links:
                    self._link_pathway_to_pathway(
                        pw_id, v["pathways"], v["direction"], v["cpd"]
                    )

    def _fix_composite_reaction_nodes(self, comp_rxn_nodes: List[Dict[str, Any]]):
        """
        Composite reaction nodes are just aggregated reactions.
        We need to change the label of these nodes from ``Pathway`` to
         ``Reaction``, and also correct their directions.

        Args:
            comp_rxn_nodes: Output of :meth:`_collect_pathways_dat_nodes`.
        """
        self.create_nodes(
            "Composite reaction",
            comp_rxn_nodes,
            self.metacyc_default_cyphers["reactions"],
        )

        # Merge the Reaction node with the Pathway node under the same ID
        for node in comp_rxn_nodes:
            self.merge_nodes("Pathway", "Reaction", "metaId", "name", node["name"])

        # Correct composite reactions
        self._set_composite_reaction_labels()

        # Fix the direction of newly-added reactions
        self._fix_reaction_direction()

    def _smiles_dat_to_graph(self, parser: MetacycParser):
        """Parse the ``atom-mappings-smiles.dat`` file.

        For details, see :meth:`setup`.
        """
        if not parser.input_files["atom_mapping"]:
            logger.warning("No atom-mappings-smiles.dat file given")
            return

        logger.info(f"Adding SMILES with {parser.input_files['atom_mapping']}")
        smiles = parser.read_smiles_dat(parser.input_files["atom_mapping"])
        all_rxns = self.get_all_nodes("Reaction", "name")
        rxn_smiles = parser.collect_atom_mapping_dat_nodes(all_rxns, smiles)
        self.add_props_to_nodes(
            "Reaction", "name", rxn_smiles, "Atom mapping SMILES", progress_bar=True
        )

    def _compounds_dat_to_graph(self, parser: MetacycParser):
        """Annotate compound nodes with the ``compounds.dat`` file.

        For details, see :meth:`setup`.
        """
        if not parser.input_files["compounds"]:
            logger.warning("No compounds.dat file given")
            return

        logger.info(f"Parsing {parser.input_files['compounds']}")
        cpd_dat = parser.read_dat_file(parser.input_files["compounds"])
        all_cpds = self.get_all_compounds()
        cpd_nodes = parser.collect_compounds_dat_nodes(all_cpds, cpd_dat)
        self.create_nodes(
            "compounds.dat",
            cpd_nodes,
            self.metacyc_default_cyphers["compounds"],
            progress_bar=True,
        )

    def _citations_dat_to_graph(self, parser: MetacycParser):
        """Annotate citation nodes with the ``pubs.dat`` file.

        For details, see :meth:`setup`.
        """
        if not parser.input_files["publications"]:
            logger.warning("No pubs.dat file given")
            return
        logger.info(
            f"Annotating publications with {parser.input_files['publications']}"
        )
        pub_dat = parser.read_dat_file(parser.input_files["publications"])

        all_cits = self.get_all_nodes("Citation", "metaId")
        cit_nodes = parser.collect_citation_dat_nodes(all_cits, pub_dat)
        self.add_props_to_nodes("Citation", "metaId", cit_nodes, "Citations")

    def _classes_dat_to_graph(self, parser: MetacycParser):
        """Parse the ``classes.dat`` file.

        Add common name and synonyms to ``Compartment`` nodes.
        Add common name, strain name, comment, and synonyms to ``Taxa`` nodes.
        """
        if not parser.input_files["classes"]:
            logger.warning("No classes.dat file given")
            return

        logger.info(f"Annotating with {parser.input_files['classes']}")
        cls_dat = parser.read_dat_file(parser.input_files["classes"])
        all_cco = self.get_all_nodes("Compartment", "name")
        if "Taxa" in self.available_node_labels:
            all_taxon = self.get_all_nodes("Taxa", "metaId")
        else:
            all_taxon = []

        cco_nodes, taxa_nodes = parser.collect_classes_dat_nodes(
            cls_dat, all_cco, all_taxon
        )
        self.add_props_to_nodes("Compartment", "name", cco_nodes, "Compartment names")
        self.add_props_to_nodes("Taxa", "metaId", taxa_nodes, "Taxa names")

    def get_all_nodes(self, label: str, prop: str) -> List[str]:
        """Fetch a property of nodes with a certain label."""
        if label not in self.available_node_labels:
            raise ValueError(f"Invalid label: {label}")
        # TODO: check for valid properties
        res = self.read(
            f"""
            MATCH (n:{label})
            RETURN DISTINCT n.{prop} AS prop;
            """
        )
        return [n["prop"] for n in res]

    def get_all_compounds(self) -> List[Tuple[str, str]]:
        """
        Fetch all ``Compound`` nodes. All Compound node with RDF have BioCyc IDs.
        """
        res = self.read(
            """
            MATCH (c:Compound)-[:hasRDF {bioQualifier: 'is'}]->(r:RDF)
            RETURN DISTINCT c.name, r.biocyc;
            """
        )  # TODO: 38 POLYMER nodes don't have BioCyc IDs
        return [(cpd["c.name"], cpd["r.biocyc"]) for cpd in res]

    def add_props_to_nodes(
        self,
        node_label: str,
        node_prop_key: str,
        nodes: List[Dict[str, Any]],
        desc: str,
        **kwargs,
    ):
        """Add properties to a list of nodes.

        Args:
            node_label: Label of the node.
            node_prop_key: Property key of the node used to locate the node.
            nodes: Properties to add to the node.
            desc: see :meth:`create_nodes`.
        """
        if node_label not in self.available_node_labels:
            logger.warning(f"Invalid label: {node_label}")
            return

        query = f"""
        UNWIND $batch_nodes AS n
          MATCH (x:{node_label} {{{node_prop_key}: n.{node_prop_key}}})
          SET x += n.props;
        """
        self.create_nodes(desc, nodes, query, **kwargs)

    def _link_pathway_to_pathway(
        self, pw: str, pws: List[str], direction: str, cpd: str
    ):
        """
        Connect a ``Pathway`` to a list of other ``Pathway``s. The ``Compound``
        that is involved in both sides of the connection is specified in the
        relationship.

        We also want to link the ``Reaction`` nodes corresponding to the
        ``Compound`` nodes that are involved in the connection. For example, if
        reaction ``R1`` produces compound ``C``, and reaction ``R2`` consumes
        compound ``C``, we want to link ``R1`` and ``R2`` with a ``isRelatedEvent``
        relationship and add the corresponding pathways as properties.
        """
        if direction == "INCOMING":
            pw_rel_type = "-[l:hasRelatedPathway]->"
            rxn_rel_type = "-[l:isRelatedEvent]->"
            l1, l2 = "hasRight", "hasLeft"
        else:
            pw_rel_type = "<-[l:hasRelatedPathway]-"
            rxn_rel_type = "<-[l:isRelatedEvent]-"
            l1, l2 = "hasLeft", "hasRight"
        self.write(
            f"""
            MATCH (pw1:Pathway {{name: $pw_id}})
            MATCH (pw2s:Pathway) WHERE pw2s.name IN $pws
            UNWIND pw2s AS pw2
            MERGE (pw1){pw_rel_type}(pw2)
                ON CREATE SET l.hasRelatedCompound = $cpd
            WITH pw1, pw2
            MATCH (pw1)-[:hasReaction]->(r1:Reaction)-[:{l1}]->(:Compound)-[:hasRDF {{bioQualifier: 'is'}}]->(:RDF {{biocyc: $cpd}})
            WITH r1, pw2
            MATCH (pw2)-[:hasReaction]->(r2:Reaction)-[:{l2}]->(:Compound)-[:hasRDF {{bioQualifier: 'is'}}]->(:RDF {{biocyc: $cpd}})
            MERGE (r1){rxn_rel_type}(r2)
                ON CREATE SET l.fromPathway = $pw_id, l.toPathway = pw2.metaId;
            """,
            pw_id=pw,
            pws=pws,
            cpd=cpd,
        )

    def _link_reaction_to_primary_compound(
        self,
        reaction_id: str,
        compound_ids: List[str],
        pathway_id: str,
        side: str,
    ):
        logger.debug(f"Reaction {reaction_id} has primary {side}s {compound_ids}")
        self.write(
            f"""
            MATCH (r:Reaction {{canonicalId: $rxn_id}})
            UNWIND $compound_ids AS cpd_id
            MATCH (cpd:Compound)-[:hasRDF {{bioQualifier: 'is'}}]->(rdf:RDF)
            WHERE rdf.biocyc = cpd_id
            MATCH (r)-[l:hasLeft|hasRight]->(cpd)
            SET l.isPrimary{side}InPathway = CASE
              WHEN l.isPrimary{side}InPathway IS NULL THEN [$pw]
              WHEN $pw IN l.isPrimary{side}InPathway THEN l.isPrimary{side}InPathway
              ELSE l.isPrimary{side}InPathway + [$pw]
              END;
            """,
            rxn_id=reaction_id,
            compound_ids=compound_ids,
            pw=pathway_id,
        )

    def merge_nodes(
        self, n1_label: str, n2_label: str, n1_attr: str, n2_attr: str, attr_val: str
    ):
        """
        Merge two nodes with the same attribute value. Note that properties is
        hard-coded to be "override", which means all node attributes of ``n2``
        will be used to override the attributes of ``n1``.

        Args:
            n1_label: The label of the first node.
            n2_label: The label of the second node.
            n1_attr: The attribute of the first node for filtering.
            n2_attr: The attribute of the second node for filtering.
            attr_val: The value of the attributes for filtering.
        """
        self.write(
            f"""
            MATCH (n1:{n1_label} {{{n1_attr}: $val}}),
                  (n2:{n2_label} {{{n2_attr}: $val}})
            CALL apoc.refactor.mergeNodes([n1, n2], {{properties: 'override'}})
            YIELD node
            RETURN node;
            """,
            val=attr_val,
        )

    def _fix_reaction_direction(self):
        """
        Seems like ``listOfReactants`` in the SBML file are always the true
        reactants in irreversible reactions, no matter if they are under
        ``LEFT`` or ``RIGHT`` in the dat file. Since the reaction nodes were all
        populated using the SBML file, we can say all reactions are either
        reversible or go from left to right.
        """
        self.write(
            """
            MATCH (r:Reaction {reactionDirection: 'right_to_left'})
            SET r.reactionDirection = 'left_to_right'
            """
        )
        self.write(
            """
            MATCH (r:Reaction {reactionDirection: 'physiol_right_to_left'})
            SET r.reactionDirection = 'physiol_left_to_right'
            """
        )

    def _delete_bad_reaction_nodes(self):
        """
        Delete reaction nodes with non-canonical metaIds and all their relationships.
        """
        self.write(
            """
            MATCH (r:Reaction) WHERE r.name <> r.canonicalId
            DETACH DELETE r;
            """
        )
        self.write(
            """
            MATCH (r:Reaction)
            REMOVE r.canonicalId;
            """
        )

    def _set_composite_reaction_labels(self):
        """
        Composite reaction nodes are labeled as ``Pathway`` and ``Reaction``.
        The ``Pathway`` label should be removed and a ``isCompositeReaction``
        property is added to the ``Reaction`` node.
        """
        self.write(
            """
            MATCH (n:Pathway:Reaction)
            REMOVE n:Pathway
            SET n.isCompositeReaction = true;
            """
        )
