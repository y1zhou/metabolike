import logging
import re
from itertools import groupby
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple, Union

import pandas as pd
from tqdm import tqdm

from metabolike.utils import add_kv_to_dict, snake_to_camel, validate_path
from .sbml import SBMLParser

logger = logging.getLogger(__name__)


class MetacycParser(SBMLParser):
    """
    Converting MetaCyc files to a Neo4j database.
    Documentation on the MetaCyc files and format FAQs can be found at:

    * MetaCyc data files download: https://metacyc.org/downloads.shtml
    * MetaCyc file formats: http://bioinformatics.ai.sri.com/ptools/flatfile-format.html
    * SBML FAQ: https://synonym.caltech.edu/documents/faq

    See :class:`.sbml.SBMLParser` for more information on SBML parsing.

    Args:
        sbml: The path to the SBML file.
        reactions: The path to the ``reaction.dat`` file. If given,
         the file will be parsed and extra annotation on ``Reaction`` nodes will
         be added.
        atom_mapping: The path to the ``atom-mappings-smiles.dat`` file. If
         given, the file will be parsed and chemical reactions in the ``SMILES``
         format will be added to the ``Reaction`` nodes.
        pathways: The path to the ``pathway.dat`` file. If given, the file will
         be parsed and pathway links will be added to the ``Reaction`` nodes.
        compounds: The path to the ``compound.dat`` file. If given, the file
         will be parsed and annotations on ``Compound`` nodes will be added.
        publications: The path to the ``publication.dat`` file. If given, the
         file will be parsed and annotations on ``Citation`` nodes will be added.
        classes: The path to the ``class.dat`` file. If given, the file will be
         parsed and annotations on ``Compartment``, ``Taxa``, and ``Compound``
         nodes will be added.

    Attributes:
        sbml_file: Filepath to the input SBML file.
        input_files: A dictionary of the paths to the input ``.dat`` files.
        missing_ids: A dictionary of sets of IDs that were not found in the
         input files. This is helpful for collecting IDs that appear to be in
         one class but are actually in another. :meth:`._report_missing_ids` can
         be used to print them out.
    """

    def __init__(
        self,
        sbml: Union[str, Path],
        reactions: Optional[Union[str, Path]] = None,
        atom_mapping: Optional[Union[str, Path]] = None,
        pathways: Optional[Union[str, Path]] = None,
        compounds: Optional[Union[str, Path]] = None,
        publications: Optional[Union[str, Path]] = None,
        classes: Optional[Union[str, Path]] = None,
    ):
        # Neo4j driver and SBML file path
        super().__init__(sbml)

        # File paths
        self.input_files = {
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

    def collect_reactions_dat_nodes(
        self, rxn_ids: Iterable[str], rxn_dat: Dict[str, List[List[str]]]
    ):
        """
        Parse entries from the reaction attribute-value file, and prepare
        nodes to add the graph database in one transaction.

        Args:
            rxn_ids: Reaction *full* ``metaId``s.
            rxn_dat: Output of :meth:`_read_dat_file` for the
             ``reaction.dat`` file.

        Returns:
            A list of dictionaries, each of which contains the information
        """
        nodes: List[Dict[str, Any]] = []
        for rxn_id in rxn_ids:
            canonical_id = self._find_rxn_canonical_id(rxn_id, rxn_dat.keys())
            if canonical_id not in rxn_dat:
                self.missing_ids["reactions"].add(canonical_id)
                continue

            node = {"name": rxn_id, "props": {"canonicalId": canonical_id}}
            lines = rxn_dat[canonical_id]
            node = self._dat_entry_to_node(
                node,
                lines,
                props_str_keys={
                    "GIBBS-0",
                    "STD-REDUCTION-POTENTIAL",
                    "REACTION-DIRECTION",
                    "REACTION-BALANCE-STATUS",
                    "SYSTEMATIC-NAME",
                    "COMMENT",  # TODO: link to other nodes
                },
                props_list_keys={"SYNONYMS", "TYPES"},
                node_list_keys={"IN-PATHWAY", "CITATIONS", "RXN-LOCATIONS"},
                prop_num_keys={"GIBBS-0", "STD-REDUCTION-POTENTIAL"},
                prop_enum_keys={"REACTION-BALANCE-STATUS", "REACTION-DIRECTION"},
            )

            nodes.append(node)

        return nodes

    def collect_pathways_dat_nodes(
        self,
        all_pws: Iterable[str],
        pw_dat: Dict[str, List[List[str]]],
        rxn_dat: Dict[str, List[List[str]]],
    ):
        pw_nodes = []
        comp_rxn_nodes = []
        parsed_pathway_ids = set()

        # We can't iterate all_pws directly because during the parsing process
        # we may add new (super-)pathways to the queue
        queue: Set[str] = set(all_pws)
        while queue:
            pw_id = queue.pop()
            # Skip if the pathway is already processed
            if pw_id in parsed_pathway_ids:
                continue
            if (pw_id not in pw_dat) and (pw_id not in rxn_dat):
                self.missing_ids["pathways"].add(pw_id)
                continue

            # pathway nodes with reaction IDs (composite reactions)
            if pw_id in rxn_dat:
                node = self.collect_reactions_dat_nodes([pw_id], rxn_dat)[0]
                # metaId should follow the SBML specification
                node["props"]["metaId"] = self._name_to_metaid(node["name"])
                if new_pws := node.get("inPathway"):
                    queue.update(new_pws)
                comp_rxn_nodes.append(node)
                parsed_pathway_ids.add(pw_id)
                continue

            # normal pathway nodes
            lines = pw_dat[pw_id]
            node = {"metaId": pw_id}
            node = self._dat_entry_to_node(
                node,
                lines,
                props_str_keys={"COMMENT", "COMMON-NAME"},
                props_list_keys={"SYNONYMS", "TYPES"},
                node_list_keys={
                    "IN-PATHWAY",
                    "SUPER-PATHWAYS",
                    "CITATIONS",
                    "SPECIES",
                    "TAXONOMIC-RANGE",
                    "PREDECESSORS",
                    "RATE-LIMITING-STEP",
                    "PRIMARY-PRODUCTS",
                    "PRIMARY-REACTANTS",
                    "REACTION-LAYOUT",
                    "PATHWAY-LINKS",
                },
            )
            if new_pws := node.get("inPathway"):
                queue.update(new_pws)
            if new_pws := node.get("superPathways"):
                queue.update(new_pws)
            pw_nodes.append(node)
            parsed_pathway_ids.add(pw_id)

        return pw_nodes, comp_rxn_nodes

    def fix_pathway_nodes(self, pw_nodes: List[Dict[str, Any]], all_rxns: Set[str]):
        """
        Some fields in the list of ``Pathway`` nodes require preprocessing
         before being fed into the database. Specifically:

        * ``predecessors`` contains a list of reaction IDs wrapped in
         parentheses. We need to extract the first ID as the target reaction,
         and take all the others as preceding events of the first one.
        * ``reactionLayout`` tells us the primary reactants and products of
         the reactions in a given pathway.
        * ``pathwayLinks`` links the pathway to other pathways through
         intermediate ``Compound``s.

        While parsing these fields, we don't add any new ``Reaction`` nodes.
        These are often hypothetical reactions and not part of the SBML file.

        Args:
            pw_nodes: Output of :meth:`_collect_pathways_dat_nodes`.
            all_rxns: All valid reaction names.
        """
        for i, n in enumerate(pw_nodes):
            if predecessors := n.get("predecessors"):
                predecessors: List[str]
                new_pred = []
                for v in predecessors:
                    if pred := self._parse_pathway_predecessors(v, all_rxns):
                        new_pred.append(pred)

                pw_nodes[i]["predecessors"] = new_pred

            if rxn_layout := n.get("reactionLayout"):
                rxn_layout: List[str]
                rxn_prim_cpds = []
                for v in rxn_layout:
                    rxn_id, d = self._parse_reaction_layout(v)
                    if rxn_id:
                        rxn_prim_cpds.append({"reaction": rxn_id, **d})

                pw_nodes[i]["reactionLayout"] = rxn_prim_cpds

            if pw_links := n.get("pathwayLinks"):
                pw_links: List[str]
                conn_pathways = []
                for v in pw_links:
                    cpd_id, pw_ids, direction = self._parse_pathway_links(v)
                    if pw_ids:
                        conn_pathways.append(
                            {"cpd": cpd_id, "pathways": pw_ids, "direction": direction}
                        )

                pw_nodes[i]["pathwayLinks"] = conn_pathways

        return pw_nodes

    def collect_atom_mapping_dat_nodes(
        self, rxn_ids: Iterable[str], smiles: Dict[str, str]
    ) -> List[Dict[str, Union[str, Dict[str, str]]]]:
        """
        Args:
            rxn_ids: The reaction name from the graph database.
            smiles: The reaction ID -> SMILES dictionary from the atom
             mapping file.
        """
        nodes = []
        for rxn_id in tqdm(rxn_ids, desc="atom_mapping.dat file"):
            canonical_id = self._find_rxn_canonical_id(rxn_id, smiles.keys())
            if canonical_id not in smiles:
                self.missing_ids["atom_mappings"].add(canonical_id)
                continue
            node = {
                "name": rxn_id,
                "props": {"smilesAtomMapping": smiles[canonical_id]},
            }
            nodes.append(node)

        return nodes

    def collect_compounds_dat_nodes(
        self, all_cpds: List[Tuple[str, str]], cpd_dat: Dict[str, List[List[str]]]
    ):
        nodes = []
        for cpd, biocyc in tqdm(all_cpds, desc="compounds.dat file"):
            # self.compounds_to_graph(cpd, biocyc, cpd_dat)
            if biocyc not in cpd_dat:
                self.missing_ids["compounds"].add(biocyc)
                continue
            lines = cpd_dat[biocyc]
            node = {"name": cpd}
            node = self._dat_entry_to_node(
                node,
                lines,
                props_str_keys={
                    "COMMON-NAME",
                    "GIBBS-0",
                    "LOGP",
                    "MOLECULAR-WEIGHT",
                    "MONOISOTOPIC-MW",
                    "POLAR-SURFACE-AREA",
                    "PKA1",
                    "PKA2",
                    "PKA3",
                    "SMILES",
                    "INCHI",
                    "COMMENT",
                    # "DBLINKS",  # TODO: parse DBLINKS
                },
                props_list_keys={"SYNONYMS"},
                node_list_keys={"CITATIONS"},
                prop_num_keys=[
                    "GIBBS-0",
                    "LOGP",
                    "MOLECULAR-WEIGHT",
                    "MONOISOTOPIC-MW",
                    "POLAR-SURFACE-AREA",
                    "PKA1",
                    "PKA2",
                    "PKA3",
                ],
            )
            nodes.append(node)

        return nodes

    def collect_citation_dat_nodes(
        self, cit_ids: Iterable[str], pub_dat: Dict[str, List[List[str]]]
    ):
        """Annotate a citation node with data from the publication.dat file.

        If there are multiple fields in the given ``cit_id``, then the fields are
        separated by colons. The first field is the citation ID, the second is the
        evidence type (in classes.dat), the third is not documented, and the
        fourth is the curator's name.

        In most cases the citation ID should match: ``PUB-[A-Z0-9]+$``,
        with a few exceptions containing double dashes, e.g. ``PUB--8``,
        or some dashes within author names, e.g. ``PUB-CHIH-CHING95``.

        Args:
            cit_ids: The citation ``metaId`` properties.
            pub_dat: The publication.dat data.
        """
        # TODO: deal with evidence frames. Evidence frames are in the form of
        # 10066805:EV-EXP-IDA:3354997111:hartmut
        nodes = []
        for cit_id in tqdm(cit_ids, desc="pubs.dat file"):
            pub_dat_id = re.sub(r"[\[\]\s,']", "", cit_id.split(":")[0].upper())
            pub_dat_id = re.sub(r"-(\d+)", r"\1", pub_dat_id)
            if pub_dat_id == "BOREJSZA-WYSOCKI94":
                pub_dat_id = "BOREJSZAWYSOCKI94"  # only exception with dash removed
            if not pub_dat_id:
                continue
            pub_dat_id = "PUB-" + pub_dat_id
            if pub_dat_id not in pub_dat:
                self.missing_ids["publications"].add(pub_dat_id)
                continue

            lines = pub_dat[pub_dat_id]
            node = {"metaId": cit_id, "props": {"citationId": pub_dat_id}}
            node = self._dat_entry_to_node(
                node,
                lines,
                props_str_keys={
                    "DOI-ID",
                    "PUBMED-ID",
                    "MEDLINE-ID",
                    "TITLE",
                    "SOURCE",
                    "YEAR",
                    "URL",
                    "REFERENT-FRAME",
                },
                props_list_keys={"AUTHORS"},
            )
            nodes.append(node)

        return nodes

    def collect_classes_dat_nodes(
        self,
        class_dat: Dict[str, List[List[str]]],
        cco_ids: Iterable[str],
        taxa_ids: Iterable[str],
    ):
        cco_nodes = []
        # Common names for cell components
        for cco in cco_ids:
            if cco not in class_dat:
                self.missing_ids["compartments"].add(cco)
                continue
            node = {"name": cco}
            node = self._dat_entry_to_node(
                node,
                class_dat[cco],
                props_str_keys={"COMMON-NAME"},
                props_list_keys={"SYNONYMS"},
            )
            cco_nodes.append(node)

        # Common names and synonyms for organisms. Some also have strain names
        taxa_nodes = []
        for taxa in tqdm(taxa_ids, desc="Taxon in classes.dat"):
            if taxa not in class_dat:
                self.missing_ids["taxon"].add(taxa)
                continue

            node = {"metaId": taxa}
            node = self._dat_entry_to_node(
                node,
                class_dat[taxa],
                props_str_keys={"COMMON-NAME", "STRAIN-NAME", "COMMENT"},
                props_list_keys={"SYNONYMS"},
                # TODO: TYPES links the taxon
            )

            taxa_nodes.append(node)

        # TODO: Evidence code in citations are in the `Evidence` attr
        return cco_nodes, taxa_nodes

    @staticmethod
    def read_dat_file(filepath: Union[str, Path]) -> Dict[str, List[List[str]]]:
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
    def read_smiles_dat(filepath: Path) -> Dict[str, str]:
        smiles_df = pd.read_table(
            filepath,
            sep="\t",
            header=None,
            names=["rxn", "smiles"],
        )
        smiles = smiles_df.set_index("rxn").to_dict()
        return smiles["smiles"]

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

        if match_canonical_id := re.match(
            r"((TRANS-)?(RXN[A-Z\d]*)-\d+)|([A-Z\d.\-+]+RXN)", rxn_id
        ):
            return match_canonical_id[0]
        else:
            raise ValueError(f"rxn_id has no canonical form: {rxn_id}")

    def _dat_entry_to_node(
        self,
        node: Dict[str, Any],
        lines: List[List[str]],
        props_str_keys: Set[str] = set(),
        props_list_keys: Set[str] = set(),
        node_str_keys: Set[str] = set(),
        node_list_keys: Set[str] = set(),
        prop_num_keys: Iterable[str] = set(),
        prop_enum_keys: Iterable[str] = set(),
    ) -> Dict[str, Any]:
        node.setdefault("props", {})
        for k, v in lines:
            if k in props_list_keys:
                add_kv_to_dict(node["props"], k, v, as_list=True)
            elif k in props_str_keys:
                add_kv_to_dict(node["props"], k, v, as_list=False)
            elif k in node_list_keys:
                add_kv_to_dict(node, k, v, as_list=True)
            elif k in node_str_keys:
                add_kv_to_dict(node, k, v, as_list=False)

        # Clean up props before writing to graph
        node["props"] = self._clean_props(
            node["props"],
            num_fields=prop_num_keys,
            enum_fields=prop_enum_keys,
        )

        return node

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
            if (k := snake_to_camel(f)) in props:
                props[k] = float(props[k])

        enum_pattern = re.compile(r"\W+")
        for f in enum_fields:
            if (k := snake_to_camel(f)) in props:
                props[k] = props[k].replace("-", "_").lower()
                props[k] = enum_pattern.sub("", props[k])

        return props

    @staticmethod
    def _name_to_metaid(s: str):
        if s[0].isnumeric():
            s = f"_{s}"
        s = s.replace("-", "__45__")
        s = s.replace(".", "__46__")
        s = s.replace("+", "")
        return s

    @staticmethod
    def _parse_pathway_predecessors(
        s: str, all_rxns: Set[str]
    ) -> Dict[str, Union[str, List[str]]]:
        """Parse the pathway predecessors.

        In most cases the string takes the form of:
         ``("RXN66-347" "1.1.1.64-RXN")``

        However, in some cases:

          * The reaction IDs are not quoted.
          * There is a single pathway ID, e.g. ``PWY-7306``.
          * There is a single reaction ID, e.g.
           ``(QUEUOSINE-TRNA-RIBOSYLTRANSFERASE-RXN)``.

        We handle the first case, and ignore the second and third cases.
        """
        if s[0] != "(":  # single pathway ID
            return {}

        if re.fullmatch(r"\(\S+\)", s):  # single reaction ID
            return {}

        s = s.replace('"', "")
        rxns = s[1:-1].split(" ")  # leading ( and trailing )
        r1, r2 = rxns[0], rxns[1:]  # could have >1 predecessors
        if r1 not in all_rxns:
            return {}
        r2 = [x for x in r2 if x in all_rxns]
        return {"r1": r1, "r2": r2}

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
            r"\(([^(]+) "  # reaction ID
            r"\(:LEFT-PRIMARIES ([^)]+)\) "  # left primary compounds
            r"\(:DIRECTION :(L2R|R2L)\) "  # reaction direction
            r"\(:RIGHT-PRIMARIES ([^)]+)\)\)"
        )
        if not (res := m.fullmatch(s)):
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
            "reactants": reactants,
            "products": products,
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
        cpd_id_rgx = re.compile(r"^\((\S+) ")
        if not (res := cpd_id_rgx.match(s)):
            logging.warning(f"Pathway links string doesn't have a compound: {s}")
            return "", [], ""
        cpd = res[1]
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
                if not (m := directed_pw_rgx.match(s)):
                    raise ValueError(f"Invalid pathway links string: {s}")
                pathways.append(m[1])
                direction = m[2]
            elif m := pw_rgx.match(s):
                pathways.append(m[0])
            else:
                raise ValueError(f"Invalid pathway links string: {s}")
            s = s[m.end() :].strip()

        # Final cleanup
        cpd = cpd.replace("|", "").replace('"', "")
        pathways = [
            p.replace('"', "") for p in pathways if (" " not in p) and p.isupper()
        ]  # remove frame IDs and quoted annotations
        direction = "INCOMING" if direction == "INCOMING" else "OUTGOING"
        return cpd, pathways, direction

    def report_missing_ids(self):
        for datfile, ids in self.missing_ids.items():
            if ids:
                logger.warning(f"The following IDs were not found in {datfile}: {ids}")
