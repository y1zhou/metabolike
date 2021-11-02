"""Functions from the brendaDb R package.

A lot of the code is translated from brendaDb:
https://bioconductor.org/packages/release/bioc/html/brendaDb.html

This implmentation focuses on extracting information from the text file, and
feeding the data into a Neo4j database.

The file is organised in an EC-number specific format. The
information on each EC-number is given in a very short and compact
way in a part of the file. The contents of each line are described
by a two/three letter acronym at the beginning of the line and
start after a TAB. Empty spaces at the beginning of a line
indicate a continuation line.

The contents are organised in ~40 information fields as given
below. Protein information is included in '#'...#', literature
citations are in '<...>', commentaries in '(...)' and field-
special information in '{...}'.

Protein information is given as the combination organism/Uniprot
accession number where available. When this information is not
given in the original paper only the organism is given.

///	indicates the end of an EC-number specific part.
"""
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Union

import pandas as pd
from lark import Lark, Token, Tree
from lark.visitors import Transformer
from met_path_finder import db

FIELDS = {
    "ACTIVATING_COMPOUND": "AC",
    "APPLICATION": "AP",
    "CLONED": "CL",
    "COFACTOR": "CF",
    "CRYSTALLIZATION": "CR",
    "ENGINEERING": "EN",
    "EXPRESSION": "EXP",
    "GENERAL_INFORMATION": "GI",
    "GENERAL_STABILITY": "GS",
    "IC50_VALUE": "IC50",
    "INHIBITORS": "IN",
    "KCAT_KM_VALUE": "KKM",  # Kcat/Km value, substrate in {...}
    "KI_VALUE": "KI",  # Ki value, inhibitor in {...}
    "KM_VALUE": "KM",  # Km value, substrate in {...}
    "LOCALIZATION": "LO",
    "METALS_IONS": "ME",
    "MOLECULAR_WEIGHT": "MW",
    "NATURAL_SUBSTRATE_PRODUCT": "NSP",  # reversibility information in {...}
    "ORGANIC_SOLVENT_STABILITY": "OSS",
    "OXIDATION_STABILITY": "OS",
    "PH_OPTIMUM": "PHO",
    "PH_RANGE": "PHR",
    "PH_STABILITY": "PHS",
    "PI_VALUE": "PI",  # isoelectric point
    "POSTTRANSLATIONAL_MODIFICATION": "PM",
    "PROTEIN": "PR",
    "PURIFICATION": "PU",
    "REACTION": "RE",
    "REACTION_TYPE": "RT",
    "RECOMMENDED_NAME": "RN",  # IUPAC accepted name
    "REFERENCE": "RF",
    "RENATURED": "REN",
    "SOURCE_TISSUE": "ST",
    "SPECIFIC_ACTIVITY": "SA",
    "STORAGE_STABILITY": "SS",
    "SUBSTRATE_PRODUCT": "SP",  # reversibility information in {...}
    "SUBUNITS": "SU",
    "SYNONYMS": "SY",
    "SYSTEMATIC_NAME": "SN",
    "TEMPERATURE_OPTIMUM": "TO",
    "TEMPERATURE_RANGE": "TR",
    "TEMPERATURE_STABILITY": "TS",
    "TURNOVER_NUMBER": "TN",  # substrate in {...}
}

FIELD_WITH_SPECIFIC_INFO = {
    "NATURAL_SUBSTRATE_PRODUCT",
    "SUBSTRATE_PRODUCT",
    "TURNOVER_NUMBER",
    "KM_VALUE",
    "KI_VALUE",
    "REFERENCE",
    "IC50_VALUE",
}

FIELD_COMMENTARY_ONLY = {
    "CLONED",
    "CRYSTALLIZATION",
    "PURIFICATION",
    "GENERAL_STABILITY",
    "STORAGE_STABILITY",
}


class Brenda:
    """Class for working with Brenda.

    Args:
        filepath: Path to the Brenda file.

    Attributes:
        db: TODO
    """

    def __init__(self):
        pass

    def parse(self, filepath: Union[str, Path]):
        # Read text file into pandas DataFrame, where the last column contains
        # the text that is to be parsed into trees.
        filepath = Path(filepath).expanduser().resolve()
        self.db = self.read_brenda(filepath)

        # TODO: parse the text into a tree

        # TODO: feed the tree into a Neo4j database

    def read_brenda(self, filepath: Path) -> pd.DataFrame:
        """Read brenda file and convert to DataFrame.

        Args:
            filepath: Path to the file

        Returns:
            A `pd.DataFrame` with columns:
                - ID: EC number, e.g. 1.1.1.1
                - field: the content of the information, e.g. protein, localization
                - description: everything else

            See :func:`clean_brenda` for details on how the data is cleaned.
        """
        lines = self._read_brenda_file(filepath)
        df = self._separate_entries(lines)
        df = self._clean_ec_number(df)

        return df

    @staticmethod
    def _read_brenda_file(filepath: Path) -> List[str]:
        """Read all non-empty lines from a file.
        Comment lines that start with '*' are ignored. The text file should be downloaded from:
        https://www.brenda-enzymes.org/download_brenda_without_registration.php

        Args:
            filepath: Path to the file

        Returns:
            List of lines in the file.
        """
        if not filepath.is_file():
            raise ValueError(f"Cannot open file: {str(filepath)}")

        res = []
        with open(filepath, "r") as f:
            for line in f:
                # Skip empty lines
                if line == "\n" or line[0] == "*":
                    continue

                # Concatenate lines that are split across multiple lines
                if line[0] == "\t":
                    res[-1] += f"{line.rstrip()}"
                else:
                    res.append(line.strip())

        return res

    @staticmethod
    def _separate_entries(lines: List[str]) -> pd.DataFrame:
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
            elif line in FIELDS:
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

    @staticmethod
    def _clean_ec_number(df: pd.DataFrame) -> pd.DataFrame:
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

    @staticmethod
    def _description_to_tree(
        text: str,
        acronym: str,
        transformer: Optional[Transformer],
        mode: str = "generic",
    ) -> Tree:
        valid_modes = {"generic", "specific_info", "reaction", "commentary"}
        if mode not in valid_modes:
            raise ValueError(f"Invalid mode: {mode}; should be one of {valid_modes}")

        # Replace placeholder with actual acronym
        grammar = fr"""
            %import .grammar.brenda_{mode} (start, _ACRONYM)
            %import .grammar.brenda_{mode} (protein_id, ref_id, description, commentary, entry)
            %import .grammar.brenda_{mode} (NUM_ID, TOKEN)
            %override _ACRONYM: "{acronym}\t"
        """

        # Parse the text into an annotated tree
        parser = Lark(grammar, start="start", parser="lalr", transformer=transformer)
        tree = parser.parse(text)

        return tree


class GenericTreeTransformer(Transformer):
    """Transform extracted values from bottom-up."""

    def NUM_ID(self, num):
        return num.update(value=int(num))

    def TOKEN(self, tok):
        return tok.update(value=tok.strip())

    def protein_id(self, children) -> Token:
        children = [x.value for x in children]
        return Token("protein_id", children)

    def ref_id(self, children) -> Token:
        children = [x.value for x in children]
        return Token("ref_id", children)

    def description(self, children) -> Token:
        res = " ".join([x.value for x in children])
        return Token("description", res)

    def commentary(self, children) -> Token:
        # protein_id, ref_id, description are the only possible cases
        res = {}
        descriptions = []
        for child in children:
            if child.type == "protein_id":
                res["protein_id"] = child.value
            elif child.type == "ref_id":
                res["ref_id"] = child.value
            else:
                descriptions.append(child.value)

        res["description"] = " ".join(descriptions)
        res["description"] = self._fix_string(res["description"])
        return Token("commentary", res)

    def entry(self, children: List[Token]) -> Dict[str, Union[str, List[int]]]:
        descriptions = []
        res = {}
        for child in children:
            if child.type == "protein_id":
                res["protein_id"] = child.value
            elif child.type == "ref_id":
                res["ref_id"] = child.value
            elif child.type == "description":
                # Concatenate description nodes
                descriptions.append(child.value)
            else:
                res["commentary"] = child.value

        res["description"] = " ".join(descriptions)
        res["description"] = self._fix_string(res["description"])
        return res

    @staticmethod
    def _fix_string(s: str):
        s = s.replace("( ", "(")
        s = s.replace(" )", ")")
        s = s.replace("[ ", "[")
        s = s.replace(" ]", "]")
        return s


if __name__ == "__main__":
    br = Brenda()
    br.parse("/mnt/data0/yi/dissertation/met_path_finder/tests/data/brenda_test.txt")
    for _, (ec_num, field, text) in br.db.iterrows():
        if field == "REACTION":
            res = br._description_to_tree(text, FIELDS[field], "reaction")
        elif field == "TRANSFERRED_DELETED":
            pass
        elif field in FIELD_WITH_SPECIFIC_INFO:
            res = br._description_to_tree(text, FIELDS[field], "specific_info")
        elif field in FIELD_COMMENTARY_ONLY:
            res = br._description_to_tree(text, FIELDS[field], "commentary")
        else:
            res = br._description_to_tree(
                text, FIELDS[field], transformer=GenericTreeTransformer()
            )
