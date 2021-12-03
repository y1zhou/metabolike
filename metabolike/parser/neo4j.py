import multiprocessing as mp
from pathlib import Path
from typing import Dict, Optional, Union

import pandas as pd
from lark import Lark
from metabolike.parser.brenda import (
    FIELDS,
    get_parser_from_field,
    read_brenda,
    text_to_tree,
)


def parse_brenda(filepath: Union[str, Path], n_jobs: int = 1, **kwargs) -> pd.DataFrame:

    """Parse the BRENDA text file into a dict.

    This implmentation focuses on extracting information from the text file,
    and feeding the data into a Neo4j database. The parser is implemented
    using Lark. A series of :class:`Transformer` classes are used to clean
    the data and convert it into the format required by Neo4j.

    Args:
        filepath: The path to the BRENDA text file.
        n_jobs: The number of parallel jobs to run.
        **kwargs: Additional keyword arguments to pass to :func:`read_brenda`.
    """
    # Read text file into pandas DataFrame, where the last column contains
    # the text that is to be parsed into trees.
    filepath = Path(filepath).expanduser().resolve()
    df = read_brenda(filepath, **kwargs)

    # Get parsers for each unique field
    parsers: Dict[str, Optional[Lark]] = {"TRANSFERRED_DELETED": None}
    for field in FIELDS.keys():
        parsers[field] = get_parser_from_field(field)

    # Parse each description into a Lark tree. To speed up parsing,
    # partition the dataframe into `n_jobs` chunks, and use dask
    # to parallelize the parsing.
    if n_jobs < 0:
        raise ValueError("n_jobs must be >= 0")

    if n_jobs == 0:
        n_jobs = mp.cpu_count()

    if n_jobs == 1:
        df["description"] = df.apply(
            lambda row: text_to_tree(row.description, parsers[row.field]),
            axis=1,
        )
    else:
        # Setup Spark
        from pyspark.sql import SparkSession

        spark = (
            SparkSession.builder.master(f"local[{n_jobs}]")
            .appName("Brenda")
            .config("spark.driver.memory", f"{n_jobs}g")
            .getOrCreate()
        )

        # Parse the text in parallel
        sdf = spark.createDataFrame(df)
        df["description"] = sdf.rdd.map(
            lambda row: text_to_tree(row.description, parsers[row.field])
        ).collect()

        spark.stop()
    return df

    # TODO: feed the tree into a Neo4j database


if __name__ == "__main__":
    parsed = parse_brenda("tests/data/brenda_test.txt", n_jobs=10, cache=True)
    print(parsed[0])
