import pandas as pd
import pandas.testing as pdt
from met_path_finder.parser.brenda import read_brenda


def test_brenda_reader():
    df = read_brenda("tests/data/brenda_test.txt", clean=True)

    assert type(df) == pd.DataFrame
    assert df.shape == (57, 3)
    pdt.assert_index_equal(df.columns, pd.Index(["ID", "field", "description"]))
