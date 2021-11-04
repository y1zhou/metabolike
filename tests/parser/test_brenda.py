from pathlib import Path

import pandas as pd
import pandas.testing as pdt
import pytest
from met_path_finder.parser.brenda import Brenda


@pytest.fixture
def brenda_data():
    br = Brenda()
    br.df = br.read_brenda(Path("tests/data/brenda_test.txt"))
    return br


def test_read_brenda(brenda_data):
    # Test that the dataframe is read correctly
    assert isinstance(brenda_data.df, pd.DataFrame)
    assert brenda_data.df.shape == (57, 3)
    pdt.assert_index_equal(
        brenda_data.df.columns, pd.Index(["ID", "field", "description"])
    )

    # DataFrame columns should match the expected values
    assert brenda_data.df["ID"].str.match(r"^(\d+\.){3}\d+$").all()
    assert brenda_data.df["field"].str.match(r"^[A-Z_05]+$").all()


def test_text_to_tree_generic(brenda_data: Brenda):
    """Test the generic grammar."""
    df = brenda_data.df
    text = (
        df.query("ID == '1.1.1.1' & field == 'PROTEIN'")
        .filter(["description"])
        .values[0, 0]
    )

    tree = brenda_data._text_to_tree(text, "PROTEIN")
    assert isinstance(tree, list)
    assert len(tree) == 164
    assert set(tree[0].keys()) - {"protein_id", "description", "ref_id"} == set()
    assert isinstance(tree[0]["protein_id"], list)
    assert isinstance(tree[0]["protein_id"][0], int)


def test_text_to_tree_reaction(brenda_data: Brenda):
    """Test the reaction grammar."""
    df = brenda_data.df
    text = (
        df.query("ID == '1.1.1.1' & field == 'SUBSTRATE_PRODUCT'")
        .filter(["description"])
        .values[0, 0]
    )

    tree = brenda_data._text_to_tree(text, "SUBSTRATE_PRODUCT")
    assert isinstance(tree, list)
    assert len(tree) == 773
    assert set(tree[0].keys()) - {"protein_id", "reaction", "ref_id"} == set()
    assert isinstance(tree[0]["protein_id"], list)
    assert isinstance(tree[0]["protein_id"][0], int)
    assert tree[0]["reaction"]["reversibility"] == "r"


def test_text_to_tree_transferred(brenda_data: Brenda):
    """Test the transferred/deleted grammar."""
    df = brenda_data.df
    text = (
        df.query("ID == '6.3.5.8' & field == 'TRANSFERRED_DELETED'")
        .filter(["description"])
        .values[0, 0]
    )
    tree = brenda_data._text_to_tree(text, "TRANSFERRED_DELETED")
    assert isinstance(tree, list)
    assert len(tree) == 1
    assert "description" in tree[0]


def test_text_to_tree_references(brenda_data: Brenda):
    """Test the references grammar."""
    df = brenda_data.df
    text = (
        df.query("ID == '1.1.1.1' & field == 'REFERENCE'")
        .filter(["description"])
        .values[0, 0]
    )

    tree = brenda_data._text_to_tree(text, "REFERENCE")
    assert isinstance(tree, list)
    assert len(tree) == 285
    assert set(tree[0].keys()) - {"ref_id", "citation", "pubmed", "paper_stat"} == set()
    assert tree[0]["pubmed"] == "6794566"
