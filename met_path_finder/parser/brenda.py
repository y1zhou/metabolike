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
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

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

FIELD_WITH_REACTIONS = {
    "REACTION",
    "NATURAL_SUBSTRATE_PRODUCT",
    "SUBSTRATE_PRODUCT",
}

FIELD_WITH_SPECIFIC_INFO = {
    "TURNOVER_NUMBER",
    "KM_VALUE",
    "KI_VALUE",
    "IC50_VALUE",
}

FIELD_COMMENTARY_ONLY = {
    "CLONED",
    "CRYSTALLIZATION",
    "PURIFICATION",
    "GENERAL_STABILITY",
    "STORAGE_STABILITY",
}

BASE_GRAMMAR = r"""
%import common.NEWLINE -> _NL

commentary : "(" _separated{content, _COMMENTARY_SEP} ")"
content    : protein_id description ref_id

!description: TOKENS+
TOKENS: TOKEN+
    | [_WS] "(" TOKEN+ ")" ["-" | ", " | _WS]  // e.g. (S)-xxx or (1 g/kg ip), parentheses but not commentary
    | [_WS] "[" /[^\]]+/  "]" // Often chemical names wrapped in [...]
protein_id : "#" _separated{NUM_ID, _LIST_SEP} "#"
ref_id     : "<" _separated{NUM_ID, _LIST_SEP} ">"

TOKEN      : [_WS] (CHAR+ | ARROW | COMPARE_NUM) [_WS]
CHAR       : /\w/
    | /[^\x00-\x7F]/  // Unicode characters
    | "?" | "." | "," | "+" | ":" | "/" | "*" | "-" | "%" | "'"
    | "(H)"  // NADP(H), NAD(H)

ARROW      : "<->" | "<-->" | "<-" | "->" | "-->"
COMPARE_NUM: /p\s?<\s?0\.05/i  // p-values
    | /[<>]/ FLOAT_NUM "%"  // percentages

_COMMENTARY_SEP : /;\s*/
_LIST_SEP  : "," | "\t"
_WS        : ("\t" | " ")+
NUM_ID     : /\d+/
FLOAT_NUM  : /\d+(\.\d+)?/

_separated{x, sep}: x (sep x)*  // A sequence of 'x sep x sep x ...'
"""


class Brenda:
    """Class for working with Brenda.

    Attributes:
        df: cleaned BRENDA text file as a pandas DataFrame.
    """

    def __init__(self):
        self.df: pd.DataFrame

    def parse(self, filepath: Union[str, Path]):
        # Read text file into pandas DataFrame, where the last column contains
        # the text that is to be parsed into trees.
        filepath = Path(filepath).expanduser().resolve()
        self.df = self.read_brenda(filepath)

        # TODO: parse the text into a tree
        for _, (ec_num, field, description) in self.df.iterrows():
            res = self._text_to_tree(description, field)

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
    def _get_parser_from_field(field: str) -> Lark:
        if field == "TRANSFERRED_DELETED":
            grammar = ""
            t: Transformer

        elif field == "REFERENCE":
            grammar = ""
            t: Transformer

        elif field in FIELD_WITH_REACTIONS:
            grammar = fr"""
                {BASE_GRAMMAR}
                %extend CHAR: "="  // Appear in reactions

                start: entry+
                entry: _ACRONYM [protein_id] reaction [commentary] [[_WS] more_commentary] [[_WS] reversibility] [[_WS] ref_id] _NL

                reaction: TOKENS+
                more_commentary:  "|" _separated{{content, _COMMENTARY_SEP}} "|"
                !reversibility:  "{{}}" | "{{r}}" | "{{ir}}"
                _ACRONYM: "{FIELDS[field]}\t"
            """
            t = ReactionTreeTransformer()

        elif field in FIELD_WITH_SPECIFIC_INFO:
            grammar = ""
            t: Transformer

        elif field in FIELD_COMMENTARY_ONLY:
            grammar = ""
            t: Transformer

        else:
            # All other cases can be parsed with the generic grammar.
            # Replace placeholder with actual acronym, and import tokens to be
            # transformed in BaseTransformer.
            grammar = fr"""
                {BASE_GRAMMAR}
                start      : entry+
                entry      : _ACRONYM [protein_id] description [commentary _WS] [ref_id] _NL
                _ACRONYM: "{FIELDS[field]}\t"
            """
            t = GenericTreeTransformer()

        return Lark(grammar, parser="lalr", transformer=t)

    def _text_to_tree(self, text: str, field: str) -> List[Dict]:
        # Get the parser for the specified mode
        parser = self._get_parser_from_field(field)

        # Parse the text into an annotated tree, then transform into dict
        tree = parser.parse(text)

        return tree.children


class BaseTransformer(Transformer):
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

    @staticmethod
    def _fix_string(s: str):
        s = re.sub(r"([\[(]) ", r"\1", s)
        s = re.sub(r" ([\])])", r"\1", s)
        s = s.replace("\t", " ")
        return s


class GenericTreeTransformer(BaseTransformer):
    """Transform extracted values from bottom-up.

    Formats the tree into a dictionary, as described in :meth:`entry`.
    """

    def description(self, children) -> Token:
        res = " ".join([x.value.strip() for x in children])
        res = self._fix_string(res)
        return Token("description", res)

    def content(self, children) -> Dict[str, Union[str, List[int]]]:
        res = {x.type: x.value for x in children}
        return res

    def commentary(self, children) -> List[Dict[str, Union[str, List[int]]]]:
        return Token("commentary", children)

    def entry(self, children: List[Token]) -> Dict[str, Any]:
        res = {x.type: x.value for x in children}
        return res


class ReactionTreeTransformer(GenericTreeTransformer):
    """
    Commentary in `(...)` are on substrates, and in `|...|` on products.
    """

    def reaction(self, children) -> Token:
        """Parse the reaction in the description node.

        There should be three parts:
            - The left-hand-side of the reaction, which is a list of chemical
                names, separated by ` + `.
            - The separator ` = `.
            - The right-hand-side of the reaction, which is a list of chemical
                names, separated by ` + `.
        """
        reaction = " ".join([x.value.strip() for x in children])
        reaction = self._fix_string(reaction)
        lhs, rhs = reaction.split(" = ")
        lhs = lhs.split(" + ")
        rhs = rhs.split(" + ")
        res = {"lhs": lhs, "rhs": rhs}
        return Token("reaction", res)

    def commentary(self, children) -> Token:
        return Token("commentary_substrate", children)

    def more_commentary(self, children) -> Token:
        return Token("commentary_product", children)

    def reversibility(self, children) -> Token:
        """
        `r` for reversible, `ir` for irreversible, `?` for unknown.
        """
        x = children[0].value
        x = x.replace("{", "")
        x = x.replace("}", "")
        if not x:
            x = "?"
        return Token("reversibility", x)

    # def entry(self, children):
    #     return children


if __name__ == "__main__":
    br = Brenda()
    br.parse("tests/data/brenda_test.txt")
