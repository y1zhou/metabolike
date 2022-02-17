import tarfile

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

