"""Functions from the brendaDb R package.

A lot of the code is translated from brendaDb:
https://bioconductor.org/packages/release/bioc/html/brendaDb.html

This implmentation focuses on extracting information from the text file, and
feeding the data into a Neo4j database.
"""
from pathlib import Path
from typing import List, Union

import pandas as pd

FIELD_SET = {
    "PROTEIN",
    "RECOMMENDED_NAME",
    "SYSTEMATIC_NAME",
    "SYNONYMS",
    "REACTION",
    "REACTION_TYPE",
    "SOURCE_TISSUE",
    "LOCALIZATION",
    "NATURAL_SUBSTRATE_PRODUCT",
    "SUBSTRATE_PRODUCT",
    "TURNOVER_NUMBER",
    "KM_VALUE",
    "PH_OPTIMUM",
    "PH_RANGE",
    "SPECIFIC_ACTIVITY",
    "TEMPERATURE_OPTIMUM",
    "TEMPERATURE_RANGE",
    "COFACTOR",
    "ACTIVATING_COMPOUND",
    "INHIBITORS",
    "KI_VALUE",
    "METALS_IONS",
    "MOLECULAR_WEIGHT",
    "POSTTRANSLATIONAL_MODIFICATION",
    "SUBUNITS",
    "PI_VALUE",
    "APPLICATION",
    "ENGINEERING",
    "CLONED",
    "CRYSTALLIZATION",
    "PURIFICATION",
    "RENATURED",
    "GENERAL_STABILITY",
    "ORGANIC_SOLVENT_STABILITY",
    "OXIDATION_STABILITY",
    "PH_STABILITY",
    "STORAGE_STABILITY",
    "TEMPERATURE_STABILITY",
    "REFERENCE",
    "IC50_VALUE",
}


def read_brenda(filepath: Union[str, Path], clean: bool = True) -> pd.DataFrame:
    """Read brenda file and convert to DataFrame.

    Args:
        filepath: Path to the file
        clean: Whether to clean the data.

    Returns:
        A `pd.DataFrame` with columns:
            - ID: EC number, e.g.
    """
    lines = read_brenda_file(filepath)
    df = separate_entries(lines)

    if clean:
        df = clean_ec_number(df)

    return df


def read_brenda_file(filepath: Union[str, Path]) -> List[str]:
    """Read all non-empty lines from a file.
    Comment lines that start with '*' are ignored. The text file should be downloaded from:
    https://www.brenda-enzymes.org/download_brenda_without_registration.php

    Args:
        filepath: Path to the file

    Returns:
        List of lines in the file.
    """
    file = Path(filepath).expanduser().resolve()
    if not file.is_file():
        raise ValueError(f"Cannot open file: {str(file)}")

    res = []
    line_need_fix = ""
    with open(file, "r") as f:
        for line in f:
            line = line.strip()
            # Skip empty lines
            if not line:
                continue

            # Concatenate lines that are split across multiple lines
            if line_need_fix:
                line = line_need_fix + line
                line_need_fix = ""

            if line[0] != "*":
                # Strip carriage returns and append to next line
                if line[-1] == "\r":
                    line_need_fix = line[:-1]
                else:
                    res.append(line)

    return res


def separate_entries(lines: List[str]) -> pd.DataFrame:
    """Convert list of lines to DataFrame.

    Args:
        lines: list of lines from :func:`read_brenda_file`.

    Returns:
        DataFrame with columns:
            - ID: EC number, e.g. 1.1.1.1
            - field: the content of the information, e.g. protein, localization
            - description: everything else
    """

    ids, fields, descriptions = [], [], []
    current_id = lines[0][3:]  # remove leading 'ID\t'
    current_field = lines[1]  # a value from FIELD_SET
    ec_info = ""

    # We skip the last line since it's just a separator '///'
    i, n = 2, len(lines) - 1
    while i < n:
        line = lines[i]

        # When we are at the end of an EC-number specific part,
        # we insert the previous entry, update ID and first field, and
        # clear the ec_info.
        if line == "///":
            ids.append(current_id)
            fields.append(current_field)
            descriptions.append(ec_info)
            current_id = lines[i + 1][3:]
            current_field = lines[i + 2]
            ec_info = ""
            i += 2

        # When we are at the end of a field, we insert the previous entry,
        # update field and clear ec_info
        elif line in FIELD_SET:
            ids.append(current_id)
            fields.append(current_field)
            descriptions.append(ec_info)

            current_field = line
            ec_info = ""

        # # Otherwise we append to the current field
        else:
            ec_info += f"{line}\n"

        i += 1

    # Insert the last entry
    ids.append(current_id)
    fields.append(current_field)
    descriptions.append(ec_info)

    return pd.DataFrame(
        {
            "ID": ids,
            "field": fields,
            "description": descriptions,
        }
    )


def clean_ec_number(df: pd.DataFrame) -> pd.DataFrame:
    """Remove deleted and transferred EC numbers.

    Some EC numbers have comments wrapped in parentheses. Most of them are
    deleted entries (in this case we remove them) or transferred entries
    (in this case we point to the new entry).

    Args:
        df Generated by :func:`read_brenda`.

    Returns:
        DataFrame with deleted and transferred entries moved to the bottom, and their:
            - `ID` being the deleted/transferred ID,
            - `field` being "TRANSFERRED_DELETED", and
            - `description` being the information included in the original ID column.
    """
    pass


if __name__ == "__main__":
    lines = read_brenda_file("~/data0/dissertation/data/brenda/brenda_2021.txt")
    df = separate_entries(lines)
    print(df)
