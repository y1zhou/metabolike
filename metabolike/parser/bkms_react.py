import tarfile
from io import StringIO
from pathlib import Path

import pandas as pd
import requests


def get_bkms_tarball(filepath: str, extract: bool = True) -> None:
    """
    The file is available at: https://bkms.brenda-enzymes.org/download.php

    The compressed file (``Reactions_BKMS.tar.gz``) includes the table in tab
    stop separated format (Excel, OpenOffice). The table contains actual data of
    BRENDA (release 2021.2, only reactions with naturally occuring substrates),
    MetaCyc (version 24.5), SABIO-RK (07/02/2021) and KEGG data, downloaded on
    the 23th of April 2012. Downloading more recent KEGG data cannot be offered
    because a KEGG license agreement would be necessary.

    Args:
        filepath: The path to store the downloaded file.
        extract: Extract the file.
    """
    url = "https://bkms.brenda-enzymes.org/download/Reactions_BKMS.tar.gz"
    r = requests.get(url, stream=True)
    with open(filepath, "wb") as f:
        for chunk in r.iter_content(chunk_size=1024):
            f.write(chunk)

    if extract:
        _extract_bkms_tarball(filepath)


def _extract_bkms_tarball(filepath: str) -> None:
    """
    Extract the tsv file from the tarball.

    Args:
        filepath: The path to the downloaded file.
    """
    tar = tarfile.open(filepath, "r:gz")
    tar.extractall()
    tar.close()


def read_bkms(filepath: str, clean: bool = True) -> pd.DataFrame:
    """
    Read the BKMS-react table and prepare it for further processing.

    The table contains random ``^M`` characters in some rows. These characters
    won't break pandas, but they will make the parsed table unexpectedly long.
    The ``clean`` parameter can be used to remove these characters.

    BRENDA takes EC numbers as identifiers, so we need entries with non-empty
    ``EC_Number`` and ``Reaction_ID_MetaCyc`` columns. There is a
    ``Reaction_ID_BRENDA`` column that's sometimes non-empty when the EC number
    is missing, but the field is not documented and it's not clear how it can be
    mapped to an entry in the BRENDA text file.

    To be extra conservative, we only keep entries with non-empty
    ``Reaction_ID_KEGG`` columns. Only reactions in MetaCyc/BioCyc that are
    associated with matching EC numbers *and* KEGG IDs will be annotated.

    Args:
        filepath: The path to the BKMS-react ``.tsv`` file.
        clean: Remove the``^M`` characters.

    Returns:
        A pandas dataframe:
    """
    fp = Path(filepath).expanduser().resolve()
    if not clean:
        df = pd.read_table(fp, sep="\t")
        return df.drop_duplicates()

    with open(fp, "r") as f:
        lines = f.readlines()
    lines = [line.replace(r"\r", "") for line in lines]
    df: pd.DataFrame = pd.read_table(StringIO("".join(lines)), sep="\t")
    df = (
        df.filter(
            [
                "EC_Number",
                "Reaction_ID_KEGG",
                "Reaction_ID_MetaCyc",
                "MetaCyc_Pathway_ID",
                "Reaction",
            ]
        )
        .query("Reaction_ID_MetaCyc.notna() & Reaction_ID_KEGG.notna()")
        .drop_duplicates()
    )
    # TODO: what are the BRENDA reaction IDs? Can we map them to MetaCyc?
    # TODO: Can the MetaCyc IDs be mapped to organism-specific BioCyc IDs?
    # TODO: clean up the EC numbers (non-canonical, missing values, etc.)
    return df
