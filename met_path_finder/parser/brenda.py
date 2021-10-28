"""Functions from the brendaDb R package.

A lot of the code is translated from brendaDb:
https://bioconductor.org/packages/release/bioc/html/brendaDb.html

This implmentation focuses on extracting information from the text file, and
feeding the data into a Neo4j database.
"""
import re
from pathlib import Path
from typing import List, Union

import pandas as pd

FIELD_SET = {
    "RECOMMENDED_NAME": "RN",
    "SYSTEMATIC_NAME": "SN",
    "REACTION": "RE",
    "PROTEIN": "PR",
    "SYNONYMS": "SY",
    "REACTION_TYPE": "RT",
    "SOURCE_TISSUE": "ST",
    "LOCALIZATION": "LO",
    "NATURAL_SUBSTRATE_PRODUCT": "NSP",
    "SUBSTRATE_PRODUCT": "SP",
    "TURNOVER_NUMBER": "TN",
    "KM_VALUE": "KM",
    "PH_OPTIMUM": "PHO",
    "PH_RANGE": "PHR",
    "SPECIFIC_ACTIVITY": "SA",
    "TEMPERATURE_OPTIMUM": "TO",
    "TEMPERATURE_RANGE": "TR",
    "COFACTOR": "CF",
    "ACTIVATING_COMPOUND": "AC",
    "INHIBITORS": "IN",
    "METALS_IONS": "ME",
    "MOLECULAR_WEIGHT": "MW",
    "POSTTRANSLATIONAL_MODIFICATION": "PM",
    "SUBUNITS": "SU",
    "PI_VALUE": "PI",
    "APPLICATION": "AP",
    "ENGINEERING": "EN",
    "CLONED": "CL",
    "CRYSTALLIZATION": "CR",
    "PURIFICATION": "PU",
    "RENATURED": "REN",
    "GENERAL_STABILITY": "GS",
    "ORGANIC_SOLVENT_STABILITY": "OSS",
    "OXIDATION_STABILITY": "OS",
    "PH_STABILITY": "PHS",
    "STORAGE_STABILITY": "SS",
    "TEMPERATURE_STABILITY": "TS",
    "REFERENCE": "RF",
    "KI_VALUE": "KI",
    "IC50_VALUE": "IC50",
}


def read_brenda(filepath: Union[str, Path], clean: bool = True) -> pd.DataFrame:
    """Read brenda file and convert to DataFrame.

    Args:
        filepath: Path to the file
        clean: Whether to clean the data.

    Returns:
        A `pd.DataFrame` with columns:
            - ID: EC number, e.g. 1.1.1.1
            - field: the content of the information, e.g. protein, localization
            - description: everything else

        See :func:`clean_brenda` for details on how the data is cleaned.
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
        See :func:`read_brenda`.
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
        df Generated by :func:`read_brenda(clean=False)`.

    Returns:
        DataFrame with deleted and transferred entries moved to the bottom, and their:
            - `ID` being the deleted/transferred ID,
            - `field` being "TRANSFERRED_DELETED", and
            - `description` being the information included in the original ID column.
    """
    df["ID"] = df["ID"].str.replace(" ()", "", regex=False)

    # Only work on entries with parentheses
    df_standard = df[~df["ID"].str.contains("(", regex=False)]
    df_nonstd = df[df["ID"].str.contains("(", regex=False)]
    df_nonstd = (
        df_nonstd.drop_duplicates("ID")
        .assign(
            field="TRANSFERRED_DELETED",
            description=lambda x: x["ID"].str.extract(r"\((.*)\)$", expand=False),
            ID=lambda x: x["ID"].str.extract(r"^((?:\d+\.){3}\d+)", expand=False),
        )
        .query("ID == ID")
    )

    return pd.concat([df_standard, df_nonstd], axis=0, ignore_index=True)


def parse_protein_num(s: str, type: str) -> str:
    """Parse protein information strings or reference strings.

    Clean a string like "#1,45, 72#" into "1,45,72".
    Consecutive commas are collapsed into one, and spaces are treated as commas.

    Args:
        s: String in the format of "#1#" or "#1,2,3#" or "<1,3>".
        type: Either "protein" or "reference".

    Returns:
        A cleaned list of protein or reference numbers separated by commas.
    """
    if not s:
        return ""
    if not type:
        raise ValueError("type must be either 'protein' or 'reference'")

    if type == "protein":
        if (not re.match(r"^#[0-9, ]+#$", s)) or re.match(r"(#,)|(,#)", s):
            raise ValueError(f"Invalid protein string: {s}")
        delim = r"#"

    elif type == "reference":
        if (not re.match(r"^<[0-9, ]+>$", s)) or re.match(r"(<,)|(,>)", s):
            raise ValueError(f"Invalid reference string: {s}")
        delim = r"[<>]"

    else:
        raise ValueError(f"Invalid type: {type}")

    s = re.sub(delim, "", s)
    s = re.sub(r"[\s,]+", ",", s)
    return s
