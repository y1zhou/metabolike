import pytest

from metabolike.utils import (
    add_kv_to_dict,
    chunk,
    generate_gene_reaction_rule,
    snake_to_camel,
    validate_path,
)


def test_validate_path():
    assert validate_path(None) is None
    f = validate_path("tests/data/pubs.dat")

    assert f.is_absolute()
    assert f.name == "pubs.dat"

    with pytest.raises(FileNotFoundError):
        validate_path("tests/data/pubs.dat.nonexistent")


def test_chunk():
    res = list(chunk([1, 2, 3, 4], 2))
    assert res == [[1, 2], [3, 4]]
    res = list(chunk([1, 2, 3, 4, 5], 2))
    assert res == [[1, 2], [3, 4], [5]]


def test_snake_to_camel():
    assert snake_to_camel("test_snake_case", "_") == "testSnakeCase"
    assert snake_to_camel("ec-code") == "ecCode"
    assert snake_to_camel("ec-code", ".") == "ec-code"
    assert snake_to_camel("kegg.reaction", ".") == "keggReaction"


def test_add_kv_to_dict():
    d = {}
    d = add_kv_to_dict(d, "k1", "v1", as_list=False)
    assert d == {"k1": "v1"}

    d = add_kv_to_dict(d, "k1", "v2", as_list=False)
    assert d == {"k1": "v2"}

    d = add_kv_to_dict(d, "k2", "v1", as_list=True)
    assert d == {"k1": "v2", "k2": ["v1"]}

    d = add_kv_to_dict(d, "k2", "v2", as_list=True)
    assert d == {"k1": "v2", "k2": ["v1", "v2"]}


def test_generate_gene_reaction_rule():
    dummy_node = {"label": "Reaction", "name": "rxn"}
    gp_nodes = [
        {
            "label": "GeneProduct",
            "metaId": f"gene__45__{i}",
            "name": f"gene-{i}",
        }
        for i in range(4)
    ]

    gp_complex = {
        "label": "GeneProductComplex",
        "metaId": "complex__45__1",
        "name": "complex-1",
    }

    gp_set = {
        "label": "GeneProductSet",
        "metaId": "set__45__1",
        "name": "set-1",
    }

    # Single gene
    rule = [(dummy_node, gp_nodes[0])]
    assert generate_gene_reaction_rule(rule) == "gene-0"

    # Single complex
    rule = [
        (dummy_node, gp_complex),
        (gp_complex, gp_nodes[0]),
        (gp_complex, gp_nodes[1]),
    ]
    assert generate_gene_reaction_rule(rule) == "gene-0 and gene-1"

    # Single set
    rule = [
        (dummy_node, gp_set),
        (gp_set, gp_nodes[0]),
        (gp_set, gp_nodes[1]),
    ]
    assert generate_gene_reaction_rule(rule) == "gene-0 or gene-1"

    # Complex within set
    rule = [
        (dummy_node, gp_set),
        (gp_set, gp_nodes[0]),
        (gp_set, gp_nodes[1]),
        (gp_set, gp_complex),
        (gp_complex, gp_nodes[2]),
        (gp_complex, gp_nodes[3]),
    ]
    assert generate_gene_reaction_rule(rule) == "gene-0 or gene-1 or ( gene-2 and gene-3 )"
