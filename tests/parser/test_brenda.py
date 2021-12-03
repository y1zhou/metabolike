from pathlib import Path

import pandas as pd
import pandas.testing as pdt
import pytest
from metabolike.parser.brenda import get_parser_from_field, read_brenda, text_to_tree


@pytest.fixture
def brenda_data():
    df = read_brenda(Path("tests/data/brenda_test.txt"))
    return df


def test_read_brenda(brenda_data):
    # Test that the dataframe is read correctly
    assert isinstance(brenda_data, pd.DataFrame)
    assert brenda_data.shape == (57, 3)
    pdt.assert_index_equal(
        brenda_data.columns, pd.Index(["ID", "field", "description"])
    )

    # DataFrame columns should match the expected values
    assert brenda_data["ID"].str.match(r"^(\d+\.){3}\d+$").all()
    assert brenda_data["field"].str.match(r"^[A-Z_05]+$").all()


def test_text_to_tree_generic(brenda_data):
    """Test the generic grammar."""
    text = (
        brenda_data.query("ID == '1.1.1.1' & field == 'PROTEIN'")
        .filter(["description"])
        .values[0, 0]
    )

    parser = get_parser_from_field("PROTEIN")
    tree = text_to_tree(text, parser)
    assert isinstance(tree, list)
    assert len(tree) == 164
    assert set(tree[0].keys()) - {"protein_id", "description", "ref_id"} == set()
    assert isinstance(tree[0]["protein_id"], list)
    assert isinstance(tree[0]["protein_id"][0], int)


def test_text_to_tree_reaction(brenda_data):
    """Test the reaction grammar."""
    text = (
        brenda_data.query("ID == '1.1.1.1' & field == 'SUBSTRATE_PRODUCT'")
        .filter(["description"])
        .values[0, 0]
    )

    parser = get_parser_from_field("SUBSTRATE_PRODUCT")
    tree = text_to_tree(text, parser)
    assert isinstance(tree, list)
    assert len(tree) == 773
    assert set(tree[0].keys()) - {"protein_id", "reaction", "ref_id"} == set()
    assert isinstance(tree[0]["protein_id"], list)
    assert isinstance(tree[0]["protein_id"][0], int)
    assert tree[0]["reaction"]["reversibility"] == "r"


def test_text_to_tree_transferred(brenda_data):
    """Test the transferred/deleted grammar."""
    text = (
        brenda_data.query("ID == '6.3.5.8' & field == 'TRANSFERRED_DELETED'")
        .filter(["description"])
        .values[0, 0]
    )

    parser = get_parser_from_field("TRANSFERRED_DELETED")
    tree = text_to_tree(text, parser)
    assert isinstance(tree, list)
    assert len(tree) == 1
    assert "description" in tree[0]


def test_text_to_tree_references(brenda_data):
    """Test the references grammar."""
    text = (
        brenda_data.query("ID == '1.1.1.1' & field == 'REFERENCE'")
        .filter(["description"])
        .values[0, 0]
    )

    parser = get_parser_from_field("REFERENCE")
    tree = text_to_tree(text, parser)
    assert isinstance(tree, list)
    assert len(tree) == 285
    assert set(tree[0].keys()) - {"ref_id", "citation", "pubmed", "paper_stat"} == set()
    assert tree[0]["pubmed"] == "6794566"


def test_text_to_tree_specific_info(brenda_data):
    """Test the specific info grammar."""
    text = (
        brenda_data.query("ID == '1.1.1.1' & field == 'TURNOVER_NUMBER'")
        .filter(["description"])
        .values[0, 0]
    )

    parser = get_parser_from_field("TURNOVER_NUMBER")
    tree = text_to_tree(text, parser)
    assert isinstance(tree, list)
    assert len(tree) == 495
    assert "substrate" in tree[0]
    assert tree[0]["description"] == "0.73"


def test_text_to_tree_commentary_only(brenda_data):
    """Test the specific info grammar."""
    text = (
        brenda_data.query("ID == '1.1.1.1' & field == 'CLONED'")
        .filter(["description"])
        .values[0, 0]
    )

    parser = get_parser_from_field("CLONED")
    tree = text_to_tree(text, parser)
    assert isinstance(tree, list)
    assert len(tree) == 85
    assert all(len(x["protein_id"]) == 1 for x in tree)
