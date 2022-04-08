from pathlib import Path
from typing import List

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
def rxn_ids(mc: MetacycParser):
    doc = mc._read_sbml(mc.sbml_file)
    model = doc.getModel()
    return [x.getMetaId() for x in model.getListOfReactions()]


def test_read_dat_file(mc: MetacycParser):
    citations = mc._read_dat_file(dat_path / "pubs.dat")

    assert isinstance(citations, dict)
    assert len(citations) == 105

    cit = citations["PUB-ARCHMICR161460"]
    assert isinstance(cit, list)
    assert all(isinstance(x, list) for x in cit)
    assert all(len(x) == 2 for x in cit)


def test_find_rxn_canonical_id(mc: MetacycParser, rxn_ids: List[str]):
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


def test_dat_entry_to_node(mc: MetacycParser):
    node = {"name": "RXN-15513", "props": {"canonicalId": "RXN-15513"}}
    rxn_dat = mc._read_dat_file(dat_path / "reactions.dat")
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
        node_list_keys={"IN-PATHWAY", "NONEXISTENT-KEY"},
        prop_num_keys={"GIBBS-0"},
        prop_enum_keys={"REACTION-BALANCE-STATUS", "REACTION-DIRECTION"},
    )

    assert isinstance(n, dict)
    assert set(n.keys()) == {"name", "props", "inPathway"}
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
    s = "(PYRUVATE  (XXX)"
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

