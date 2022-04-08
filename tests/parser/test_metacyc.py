from pathlib import Path
from typing import Dict, List

import pytest
from metabolike.parser import MetacycParser

dat_path = Path("tests/data/")


@pytest.fixture(name="mc")
def mc_parser():
    # No need to pass an actual database to the parser
    return MetacycParser(
        None,
        sbml=dat_path / "metabolic-reactions.sbml",
        reactions=dat_path / "reactions.dat",
        atom_mapping=dat_path / "atom-mappings-smiles.dat",
        pathways=dat_path / "pathways.dat",
        compounds=dat_path / "compounds.dat",
        publications=dat_path / "pubs.dat",
        classes=dat_path / "classes.dat",
    )


@pytest.fixture
def rxn_dat(mc: MetacycParser):
    return mc._read_dat_file(mc.input_files["reactions"])


def test_read_dat_file(mc: MetacycParser):
    citations = mc._read_dat_file(mc.input_files["publications"])

    assert isinstance(citations, dict)
    assert len(citations) == 105

    cit = citations["PUB-ARCHMICR161460"]
    assert isinstance(cit, list)
    assert all(isinstance(x, list) for x in cit)
    assert all(len(x) == 2 for x in cit)


def test_read_smiles_dat(mc: MetacycParser):
    smiles = mc._read_smiles_dat(mc.input_files["atom_mapping"])

    assert isinstance(smiles, dict)
    assert len(smiles) == 12


def test_find_rxn_canonical_id(mc: MetacycParser):
    doc = mc._read_sbml(mc.sbml_file)
    model = doc.getModel()
    rxn_ids = [x.getMetaId() for x in model.getListOfReactions()]

    # Valid reaction IDs
    assert mc._find_rxn_canonical_id("RXN-15513", rxn_ids) == "RXN-15513"
    assert mc._find_rxn_canonical_id("6PFRUCTPHOS-RXN", rxn_ids) == "6PFRUCTPHOS-RXN"
    assert (
        mc._find_rxn_canonical_id(
            "F16BDEPHOS-RXN[CCO-CYTOSOL]-FRUCTOSE-16-DIPHOSPHATE/WATER//FRUCTOSE-6P/Pi.59.",
            rxn_ids,
        )
        == "F16BDEPHOS-RXN"
    )

    # Invalid reaction IDs
    with pytest.raises(ValueError):
        mc._find_rxn_canonical_id("RXN0", rxn_ids)

    # Return whatever if found in rxn_ids
    assert mc._find_rxn_canonical_id("RXN0", ["RXN0"] + rxn_ids) == "RXN0"


def test_clean_props(mc: MetacycParser):
    props = {
        "name": "RXN-15513",
        "gibbs0": "-1.3200073",
        "reactionBalanceStatus": ":BALANCED",
    }

    assert mc._clean_props(
        props, num_fields=["GIBBS-0"], enum_fields=["REACTION-BALANCE-STATUS"]
    ) == {
        "name": "RXN-15513",
        "gibbs0": -1.3200073,
        "reactionBalanceStatus": "balanced",
    }


def test_dat_entry_to_node(mc: MetacycParser, rxn_dat: Dict[str, List[List[str]]]):
    node = {"name": "RXN-15513", "props": {"canonicalId": "RXN-15513"}}
    lines = rxn_dat["RXN-15513"]

    n = mc._dat_entry_to_node(node, lines)
    assert n == node

    n = mc._dat_entry_to_node(
        node,
        lines,
        props_str_keys={
            "GIBBS-0",
            "REACTION-DIRECTION",
            "REACTION-BALANCE-STATUS",
        },
        props_list_keys={"TYPES"},
        node_str_keys={"EC-NUMBER"},
        node_list_keys={"IN-PATHWAY", "NONEXISTENT-KEY"},
        prop_num_keys={"GIBBS-0"},
        prop_enum_keys={"REACTION-BALANCE-STATUS", "REACTION-DIRECTION"},
    )

    assert isinstance(n, dict)
    assert set(n.keys()) == {"name", "props", "inPathway", "ecNumber"}
    assert set(n["props"].keys()) == {
        "canonicalId",
        "gibbs0",
        "reactionDirection",
        "reactionBalanceStatus",
        "types",
    }

    assert isinstance(n["inPathway"], list)
    assert len(n["inPathway"]) == 7
    assert isinstance(n["props"]["gibbs0"], float)
    assert n["props"]["reactionDirection"] == "reversible"


def test_name_to_metaid(mc: MetacycParser):
    assert mc._name_to_metaid("RXN-15513") == "RXN__45__15513"
    assert mc._name_to_metaid("1.1.1.1-RXN") == "_1__46__1__46__1__46__1__45__RXN"
    assert mc._name_to_metaid("FE2+-RXN") == "FE2__45__RXN"


def test_parse_pathway_predecessors(mc: MetacycParser):
    all_rxns = {"RXN-13760", "RXN-19880", "PGLUCISOM-RXN"}

    # single pathway ID
    assert mc._parse_pathway_predecessors("PWY-1", all_rxns) == {}

    # single reaction ID
    assert mc._parse_pathway_predecessors("(RXN-13760)", all_rxns) == {}
    assert mc._parse_pathway_predecessors("(RXN-1)", all_rxns) == {}

    # Normal cases
    pred = mc._parse_pathway_predecessors("(RXN-19880 RXN-13760)", all_rxns)
    assert pred == {"r1": "RXN-19880", "r2": ["RXN-13760"]}
    pred = mc._parse_pathway_predecessors('("RXN-19880" "RXN-13760")', all_rxns)
    assert pred == {"r1": "RXN-19880", "r2": ["RXN-13760"]}
    s = "(RXN-19880 RXN-13760 PGLUCISOM-RXN)"
    assert mc._parse_pathway_predecessors(s, all_rxns) == {
        "r1": "RXN-19880",
        "r2": ["RXN-13760", "PGLUCISOM-RXN"],
    }

    # r1 not in all_rxns
    pred = mc._parse_pathway_predecessors("(RXN-1 RXN-13760)", all_rxns)
    assert pred == {}

    # r2 not in all_rxns
    pred = mc._parse_pathway_predecessors("(RXN-19880 RXN-13760 RXN-1)", all_rxns)
    assert pred == {"r1": "RXN-19880", "r2": ["RXN-13760"]}


def test_parse_reaction_layout(mc: MetacycParser):
    s = (
        "(PGLUCISOM-RXN "
        "(:LEFT-PRIMARIES D-glucopyranose-6-phosphate) "
        "(:DIRECTION :L2R) "
        "(:RIGHT-PRIMARIES FRUCTOSE-6P))"
    )

    assert mc._parse_reaction_layout(f"{s} ") == ("", {})
    rxn_id, d = mc._parse_reaction_layout(s)
    assert rxn_id == "PGLUCISOM-RXN"
    assert d == {
        "reactants": ["D-glucopyranose-6-phosphate"],
        "products": ["FRUCTOSE-6P"],
    }
    s = s.replace(":L2R", ":R2L")
    _, d = mc._parse_reaction_layout(s)
    assert d == {
        "reactants": ["FRUCTOSE-6P"],
        "products": ["D-glucopyranose-6-phosphate"],
    }


def test_parse_pathway_links(mc: MetacycParser):
    # invalid string
    s = "PYRUVATE ALANINE-SYN2-PWY"
    assert mc._parse_pathway_links(s) == ("", [], "")

    s = "(PYRUVATE  ))"
    with pytest.raises(ValueError):
        mc._parse_pathway_links(s)
    s = "(PYRUVATE (XXX)"
    with pytest.raises(ValueError):
        mc._parse_pathway_links(s)

    # simplest case
    s = "(PYRUVATE ALANINE-SYN2-PWY)"
    cpd, pathways, direction = mc._parse_pathway_links(s)
    assert cpd == "PYRUVATE"
    assert pathways == ["ALANINE-SYN2-PWY"]
    assert direction == "OUTGOING"

    # multiple pathways
    s = "(PYRUVATE PWY-5481 PWY-7351)"
    cpd, pathways, direction = mc._parse_pathway_links(s)
    assert cpd == "PYRUVATE"
    assert pathways == ["PWY-5481", "PWY-7351"]
    assert direction == "OUTGOING"

    # directed link
    s = "(GLYCEROL (TRANS-RXN-131 . :INCOMING))"
    cpd, pathways, direction = mc._parse_pathway_links(s)
    assert cpd == "GLYCEROL"
    assert pathways == ["TRANS-RXN-131"]
    assert direction == "INCOMING"

    s = '("D-glucopyranose" PWY-6724 (RXN-14352 . :INCOMING))'
    cpd, pathways, direction = mc._parse_pathway_links(s)
    assert cpd == "D-glucopyranose"
    assert pathways == ["PWY-6724", "RXN-14352"]
    assert direction == "INCOMING"

    # frame ID
    s = "(|Starch| PWY-6724)"
    cpd, pathways, direction = mc._parse_pathway_links(s)
    assert cpd == "Starch"
    assert pathways == ["PWY-6724"]
    assert direction == "OUTGOING"

    s = "(G3P |Biosynthesis|)"
    cpd, pathways, direction = mc._parse_pathway_links(s)
    assert cpd == "G3P"
    assert pathways == []
    assert direction == "OUTGOING"

    # random string
    s = '(ACETYL-P "conversion to acetate")'
    cpd, pathways, direction = mc._parse_pathway_links(s)
    assert cpd == "ACETYL-P"
    assert pathways == []
    assert direction == "OUTGOING"


def test_collect_reactions_dat_nodes(
    mc: MetacycParser, rxn_dat: Dict[str, List[List[str]]]
):
    rxn_ids = ["RXN-15513", "F16BDEPHOS-RXN"]
    rxn_nodes = mc._collect_reactions_dat_nodes(rxn_ids, rxn_dat)
    assert isinstance(rxn_nodes, list)
    assert len(rxn_nodes) == 2
    assert len(mc.missing_ids["reactions"]) == 0

    rxn_ids = ["RXN-FAKE"]
    with pytest.raises(ValueError):
        mc._collect_reactions_dat_nodes(rxn_ids, rxn_dat)

    rxn_ids = ["RXN-1"]
    _ = mc._collect_reactions_dat_nodes(rxn_ids, rxn_dat)
    assert len(mc.missing_ids["reactions"]) == 1


def test_collect_pathways_dat_nodes(
    mc: MetacycParser, rxn_dat: Dict[str, List[List[str]]]
):
    pw_dat = mc._read_dat_file(mc.input_files["pathways"])
    pws = ["PWY-3801", "PWY-3801"]
    pw_nodes, comp_rxn_nodes = mc._collect_pathways_dat_nodes(pws, pw_dat, rxn_dat)

    assert isinstance(pw_nodes, list)
    assert [x["metaId"] for x in pw_nodes] == ["PWY-3801", "PWY-7345"]
    assert isinstance(comp_rxn_nodes, list)
    assert len(comp_rxn_nodes) == 0
    assert len(mc.missing_ids["pathways"]) == 0

    # test with missing pathway ID
    _ = mc._collect_pathways_dat_nodes(["PWY-FAKE"], pw_dat, rxn_dat)
    assert len(mc.missing_ids["pathways"]) == 1

    # test with fake composite reaction node
    pws.append("RXN-15513")
    pw_nodes, comp_rxn_nodes = mc._collect_pathways_dat_nodes(pws, pw_dat, rxn_dat)
    assert len(comp_rxn_nodes) == 1
    assert len(pw_nodes) == 15


def test_fix_pathway_nodes(mc: MetacycParser):
    n = {
        "metaId": "UDPNAGSYN-PWY",
        "predecessors": [
            '("2.3.1.157-RXN" "5.4.2.10-RXN")',
            '("NAG1P-URIDYLTRANS-RXN" "2.3.1.157-RXN")',
        ],
        "reactionLayout": [
            "(2.3.1.157-RXN (:LEFT-PRIMARIES GLUCOSAMINE-1P) (:DIRECTION :L2R) (:RIGHT-PRIMARIES N-ACETYL-D-GLUCOSAMINE-1-P))",
            "(5.4.2.10-RXN (:LEFT-PRIMARIES D-GLUCOSAMINE-6-P) (:DIRECTION :L2R) (:RIGHT-PRIMARIES GLUCOSAMINE-1P))",
            "(NAG1P-URIDYLTRANS-RXN (:LEFT-PRIMARIES N-ACETYL-D-GLUCOSAMINE-1-P) (:DIRECTION :L2R) (:RIGHT-PRIMARIES UDP-N-ACETYL-D-GLUCOSAMINE))",
        ],
        "pathwayLinks": ["(FRUCTOSE-6P PWY-5484)"],
    }

    all_rxns = {"2.3.1.157-RXN", "5.4.2.10-RXN"}
    mc._fix_pathway_nodes([n], all_rxns)  # dict is mutable

    assert isinstance(n["predecessors"], list)
    assert isinstance(n["reactionLayout"], list)
    assert isinstance(n["pathwayLinks"], list)

    assert len(n["predecessors"]) == 1  # dropped NAG1P-URIDYLTRANS-RXN
    assert len(n["reactionLayout"]) == 3
    assert len(n["pathwayLinks"]) == 1


def test_collect_atom_mapping_dat_nodes(mc: MetacycParser):
    smiles = mc._read_smiles_dat(mc.input_files["atom_mapping"])

    nodes = mc._collect_atom_mapping_dat_nodes({"RXN-15513", "PEPSYNTH-RXN"}, smiles)

    assert isinstance(nodes, list)
    assert len(nodes) == 2
    assert "props" in nodes[0]
    assert "smilesAtomMapping" in nodes[0]["props"]

    _ = mc._collect_atom_mapping_dat_nodes({"RXN-000"}, smiles)
    assert mc.missing_ids["atom_mappings"] == {"RXN-000"}


def test_collect_compounds_dat_nodes(mc: MetacycParser):
    cpds_dat = mc._read_dat_file(mc.input_files["compounds"])
    cpds = [
        ("phosphate", "Pi"),
        ("pyruvate", "PYRUVATE"),
        ("2-phospho-D-glycerate", "2-PG"),
    ]
    nodes = mc._collect_compounds_dat_nodes(cpds, cpds_dat)
    assert isinstance(nodes, list)
    assert len(nodes) == 3
    assert len(mc.missing_ids["compounds"]) == 0

    cpds = [("fake", "fake")]
    _ = mc._collect_compounds_dat_nodes(cpds, cpds_dat)
    assert len(mc.missing_ids["compounds"]) == 1


def test_collect_citation_dat_nodes(mc: MetacycParser):
    cit_dat = mc._read_dat_file(mc.input_files["publications"])

    test_ids = [
        "ARCHMICR161460",
        "16672461",
        "BOREJSZA-WYSOCKI94",
        "CHIH-CHING95",
        "[ ]",
    ]
    nodes = mc._collect_citation_dat_nodes(test_ids, cit_dat)

    assert isinstance(nodes, list)
    assert len(nodes) == 2
    assert mc.missing_ids["publications"] == {
        "PUB-BOREJSZAWYSOCKI94",
        "PUB-CHIH-CHING95",
    }


def test_collect_classes_dat_nodes(mc: MetacycParser):
    classes_dat = mc._read_dat_file(mc.input_files["classes"])
    cpt_ids = ["CCO-CHROM-STR", "CCO-CYTOSOL", "CCO-IN"]
    taxa_ids = ["TAX-9606", "TAX-10090", "TAX-0"]

    cpt_nodes, taxa_nodes = mc._collect_classes_dat_nodes(
        classes_dat, cpt_ids, taxa_ids
    )

    assert mc.missing_ids["compartments"] == {"CCO-IN"}
    assert mc.missing_ids["taxon"] == {"TAX-0"}

    assert isinstance(cpt_nodes, list)
    assert len(cpt_nodes) == 2
    assert isinstance(taxa_nodes, list)
    assert len(taxa_nodes) == 2
