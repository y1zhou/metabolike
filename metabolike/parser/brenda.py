"""Reading and parsing the BRENDA text file.

A lot of the code is translated from the `brendaDb R package`_.
The `BRENDA text file`_ is organised in an EC-number specific format. The
information on each EC-number is given in a very short and compact
way in a part of the file. The contents of each line are described
by a two/three letter acronym at the beginning of the line and
start after a TAB. Empty spaces at the beginning of a line
indicate a continuation line.

.. _`brendaDb R package`: https://bioconductor.org/packages/release/bioc/html/brendaDb.html
.. _`BRENDA text file`: https://www.brenda-enzymes.org/download_brenda_without_registration.php

The contents are organised in ~40 information fields as given
below. Protein information is included in ``#...#``, literature
citations are in ``<...>``, commentaries in ``(...)`` and field-
special information in ``{...}``.

It's not officially documented, but some fields also have commentaries
wrapped in ``|...|``. These are usually for reaction-related fields, so
``(...)`` would be for substrates and ``|...|`` for products.

Protein information is given as the combination organism/Uniprot
accession number where available. When this information is not
given in the original paper only the organism is given.

``///``	indicates the end of an EC-number specific part.
"""

from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd
from lark import Lark
from metabolike.parser.brenda_transformer import *

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

# The following fields have additional commentary information,
# reversibility, and a specific reaction format in the descriptions.
FIELD_WITH_REACTIONS = {
    "REACTION",
    "NATURAL_SUBSTRATE_PRODUCT",
    "SUBSTRATE_PRODUCT",
}

# The following fields have substrate or inhibitor information
# wrapped in `{...}`.
FIELD_WITH_SPECIFIC_INFO = {
    "TURNOVER_NUMBER",
    "KCAT_KM_VALUE",
    "KM_VALUE",
    "KI_VALUE",
    "IC50_VALUE",
}

# The following fields don't have descriptions; there is only
# a blob of text wrapped in `(...)`.
FIELD_COMMENTARY_ONLY = {
    "CLONED",
    "CRYSTALLIZATION",
    "PURIFICATION",
    "GENERAL_STABILITY",
    "RENATURED",
    "STORAGE_STABILITY",
}

BASE_GRAMMAR = r"""
%import common.NEWLINE -> _NL

commentary : /\((?=#)/ _separated{content, _COMMENTARY_SEP} /(?<=>)\)/
content    : protein_id description ref_id

!description : TOKENS+
protein_id   : "#" _separated{NUM_ID, _LIST_SEP} "#"
ref_id       : "<" _separated{NUM_ID, _LIST_SEP} ">"

TOKENS: TOKEN+
TOKEN : [_WS] (CHAR+ | INLINE_CHEMICAL) [_WS]
CHAR  : /\w/
      | /[^\x00-\x7F]/  // Unicode characters
      | ARROW | COMPARE_NUM
      | "?" | "." | "," | "+" | ":" | "-" | "/" | "!"  // Common punctuation
      | "*" | "%" | "'" | "\"" | "&" | ";" | "~" | "`" | "\\"  // Rare symbols
      | "[" | "]"  // Wraps chemical names
      | "="  // Appears in REACTION and very rarely in PROTEIN

ARROW          : "<->" | "<-->" | "<-" | "->" | "-->" | "<<" | ">" | /<\t?>/ | /\s?<\s/
               | /<(?=[A-Za-z])/  // Comparison of chemicals in commentary
COMPARE_NUM    : /p\s?<\s?0\.05/i  // p-values

//INLINE_CHEMICAL: /\(+(?!#)/ (CHAR | _WS | COMPARE_NUM | ARROW | "(" | ")")+ ")"
INLINE_CHEMICAL: /\((?!#\d)/ | /(?<!\d>)\)/

_COMMENTARY_SEP : /;\s*(?=#)/
_LIST_SEP       : "," | "\t"
_WS             : ("\t" | " ")+
NUM_ID          : /\d+/
FLOAT_NUM       : NUM_ID "." NUM_ID

_separated{x, sep}: x (sep x)*  // A sequence of 'x sep x sep x ...'

start: entry+
"""


def get_parser_from_field(field: str) -> Optional[Lark]:
    if field == "TRANSFERRED_DELETED":
        return None

    elif field == "REFERENCE":
        grammar = fr"""
            {BASE_GRAMMAR}
            entry      : _ACRONYM ref_id citation [_WS] pubmed [_WS paper_stat] _NL
            citation   : /.+(?={{Pubmed:)/
            pubmed     : /\{{Pubmed:\d*\}}+/
            paper_stat : "(" _separated{{STATS, ","}} ")"
            STATS      : /[\w\/.\-\s]+/
            _ACRONYM: "{FIELDS[field]}\t"
        """
        t = RefTreeTransformer()

    elif field in FIELD_WITH_REACTIONS:
        grammar = fr"""
            {BASE_GRAMMAR}
            entry: _ACRONYM [protein_id] reaction [commentary] [[_WS] more_commentary] [[_WS] reversibility] [[_WS] ref_id] _NL

            reaction: TOKENS+
            %extend CHAR: /\{{(?!(r|ir|\?)?\}})/ | "}}"  // inline curly brackets around chemicals
            %override content: protein_id description [[_WS] reversibility [_WS]] ref_id  // reversibility in commentary

            more_commentary:  "|" _separated{{content, _COMMENTARY_SEP}} "|"
            !reversibility:  "{{}}" | "{{r}}" | "{{ir}}" | "{{?}}"
            _ACRONYM: "{FIELDS[field]}\t"
        """
        t = ReactionTreeTransformer()

    elif field in FIELD_WITH_SPECIFIC_INFO:
        grammar = fr"""
            {BASE_GRAMMAR}
            entry      : _ACRONYM protein_id description substrate [_WS commentary] [_WS] ref_id _NL

            %extend CHAR: /(?<=[^\s])[\{{\}}](?=[^\s])/
            substrate: [_WS] "{{" (TOKENS | "{{" | "}}" | "(" | ")")+ "}}"
            _ACRONYM: "{FIELDS[field]}\t"
        """
        t = SpecificInfoTreeTransformer()

    elif field in FIELD_COMMENTARY_ONLY:
        grammar = fr"""
            {BASE_GRAMMAR}
            entry      : _ACRONYM protein_id _WS [description] ref_id _NL
            %override description: /\(.*\)/ _WS
            _ACRONYM: "{FIELDS[field]}\t"
        """
        t = CommentaryOnlyTreeTransformer()

    else:
        # All other cases can be parsed with the generic grammar.
        # Replace placeholder with actual acronym, and import tokens to be
        # transformed in BaseTransformer.
        grammar = fr"""
            {BASE_GRAMMAR}
            %extend CHAR: "{{" | "}}" // Appears in ENGINEERING and INHIBITORS

            entry      : _ACRONYM [protein_id] description [commentary] [_WS] [ref_id] _NL
            _ACRONYM: "{FIELDS[field]}\t"
        """
        t = GenericTreeTransformer()

    return Lark(grammar, parser="lalr", transformer=t)


def text_to_tree(text: str, parser: Optional[Lark]) -> List[Dict]:
    if parser is None:
        # Simply return the text for TRANSFERRED_DELETED fields
        return [{"description": text}]

    # Parse the text into an annotated tree, then transform into dict
    tree = parser.parse(text)
    return tree.children


def read_brenda(filepath: Path, cache: bool = False) -> pd.DataFrame:
    """Read brenda file and convert to DataFrame.

    Args:
        filepath: Path to the file

    Returns:
        A :class:`pandas.DataFrame` with columns:

           * ID: EC number, e.g. 1.1.1.1
           * field: the content of the information, e.g. protein, localization
           * description: everything else
    """
    cache_file = filepath.with_suffix(".parquet")
    if cache_file.exists() and cache:
        df = pd.read_parquet(cache_file)
        return df

    lines = _read_brenda_file(filepath)
    df = _separate_entries(lines)
    df = _clean_ec_number(df)

    # Fix spefical cases that break the parser
    regexes = [
        ("<Swissprot>", "", False),  # SYNONYMS has this in ref_id
        ("&#8242;", "'", False),  # &#8242; is a single quote
        ("Rmyt2<", "Rmyt2", False),  # Extra < in string
        (
            r"[\x00-\x08\x0b-\x1f\x7f-\x9f]",
            " ",
            True,
        ),  # Remove all control characters other than \t and \n
        # Remove # in text that's not protein ID
        (r"(?<=\w)#(?=\w)", "", True),
        (r"\s#([A-Z])", r" \1", True),
        (r"#(\d+[,\t]?[^\d,\t#])", r"number \1", True),
        (r"([^\s\d\|\(])#", r"\1", True),
        (r"#2#\s#> <38>;", "", True),  # 3.4.21.69
        (
            r"(stereoselectivity) :(.*<5,6,7>) (<3,6>;)",
            r"\1 \3 \2; ",
            True,
        ),  # nested commentary in 1.1.1.270
        (
            r">(4#.+<20>) (<29>)",
            r"\2; #\1",
            True,
        ),  # nested commentary in 2.7.9.1 REACTION
        ("<6> <17>", "<6,17>", False),  # nested commentary in 3.4.21.79 NSP
        ("#2# <4>;", "<4>;", False),  # nested commentary in 5.1.3.3 REACTION
        # Missing opening/closing parenthees
        (r"\sE\)-farnesyl", " (E)-farnesyl", True),
        ("(4a-", "4a-", False),  # (4a-hydroxytetrahydrobiopterin
        ("SN\t2E,6E)", "SN\t(2E,6E)", False),
        (" dU)230 ", " (dU)230 ", False),  # 4.2.99.18 SP
        (
            "carboxyethylsulfanylthiocarbonylamino\t",
            "carboxyethylsulfanylthiocarbonylamino)",
            False,
        ),  # 1.8.1.7 INHIBITORS,
        ("1->4>)", "1->4)", False),  # 4.2.2.8 SUBSTRATE_PRODUCT
        # Replace < with "less than" for non-ref ID cases
        (r"<(\d+[^,\d>])", r"less than \1", True),
        (r"<(\d+,\s\d+)", r"less than \1", True),
        (
            "#1,2,4-8,12-14#",
            "",
            False,
        ),  # Ranges in protein ID for 4.2.1.49
    ]
    for x in regexes:
        df.description = df.description.str.replace(x[0], x[1], regex=x[2])

    df = df[df.description != "SN\n"]  # empty systamatic name

    # Remove redundant `|` in description
    mask = df.description.str.contains(r"\|[^#]") & (
        ~df.field.isin(FIELD_WITH_REACTIONS)
    )
    df.loc[mask, "description"] = df.loc[mask, "description"].str.replace(
        "|", "", regex=False
    )

    # Irregular commentary in 2.2.1.3 REACTION
    mask = (df.ID == "2.2.1.3") & (df.field == "REACTION")
    df.loc[mask, "description"] = df.loc[mask, "description"].str.replace(
        r"\(#.*$", "", regex=True
    )

    df.reset_index(drop=True, inplace=True)

    if cache:
        df.to_parquet(str(cache_file), index=False)

    return df


def parse_brenda(filepath: Union[str, Path], **kwargs) -> pd.DataFrame:

    """Parse the BRENDA text file into a dict.

    This implmentation focuses on extracting information from the text file,
    and feeding the data into a Neo4j database. The parser is implemented
    using Lark. A series of :class:`lark.visitors.Transformer` classes are used
    to clean the data and convert it into the format required by Neo4j.

    Args:
        filepath: The path to the BRENDA text file.
        **kwargs: Additional keyword arguments to pass to :func:`read_brenda`.

    Returns:
        A :class:`pandas.DataFrame` with the ``description`` column from
        :func:`read_brenda` transformed into lists of dicts.
    """
    # Read text file into pandas DataFrame, where the last column contains
    # the text that is to be parsed into trees.
    filepath = Path(filepath).expanduser().resolve()
    df = read_brenda(filepath, **kwargs)

    # Get parsers for each unique field
    parsers: Dict[str, Optional[Lark]] = {"TRANSFERRED_DELETED": None}
    for field in FIELDS.keys():
        parsers[field] = get_parser_from_field(field)

    df["description"] = df.apply(
        lambda row: text_to_tree(row.description, parsers[row.field]),
        axis=1,
    )
    return df

    # TODO: feed the tree into a Neo4j database


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


def _separate_entries(lines: List[str]) -> pd.DataFrame:
    """Convert list of lines to DataFrame.

    Args:
        lines: list of lines from :func:`_read_brenda_file`.

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


def _clean_ec_number(df: pd.DataFrame) -> pd.DataFrame:
    """Remove deleted and transferred EC numbers.

    Some EC numbers have comments wrapped in parentheses. Most of them are
    deleted entries (in this case we remove them) or transferred entries
    (in this case we point to the new entry).

    Args:
        df Generated by :func:`read_brenda`.

    Returns:
        :class:`pandas.DataFrame` with deleted and transferred entries moved to
        the bottom, and their:

           * ``ID`` being the deleted/transferred ID,
           * ``field`` being "TRANSFERRED_DELETED", and
           * ``description`` being the information included in the original ID column.
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
