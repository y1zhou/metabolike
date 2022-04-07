import logging
from itertools import islice
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple, Union, overload

logger = logging.getLogger(__name__)


@overload
def validate_path(filepath: Union[str, Path]) -> Path:
    ...


@overload
def validate_path(filepath: None) -> None:
    ...


def validate_path(filepath: Optional[Union[str, Path]]) -> Optional[Path]:
    if not filepath:
        return None
    f = Path(filepath).expanduser().resolve()
    if not f.is_file():
        logger.error(f"File does not exist: {f}")
        raise FileNotFoundError(str(f))

    return f


def chunk(itr: Iterable, chunk_size: int) -> Iterable:
    """Chunk an iterator into chunks of size chunk_size.

    Ref: https://stackoverflow.com/a/71572942/5925357
    """
    itr = iter(itr)
    while slice := list(islice(itr, chunk_size)):
        yield slice


def snake_to_camel(s: str, sep: str = "-") -> str:
    """Convert snake_case to CamelCase.

    Args:
        s: String to convert.
        sep: Separator to use.

    Returns:
        camelCase string.
    """
    s_camel = [x.capitalize() for x in s.split(sep)]
    s_camel[0] = s_camel[0].lower()
    return "".join(s_camel)


def add_kv_to_dict(
    d: Dict[str, Any], k: str, v: Any, as_list: bool = False
) -> Dict[str, Any]:
    """Add a key-value pair to a dictionary.

    Args:
        d: Dictionary to add the key-value pair to.
        k: Key to add.
        v: Value to add.

    Returns:
        Dictionary with the key-value pair added. If the key already exists,
        the value is saved as a list.
    """
    k_camel = snake_to_camel(k)
    if as_list:
        d.setdefault(k_camel, []).append(v)
    else:
        d[k_camel] = v

    return d


def generate_gene_reaction_rule(genes: List[Tuple[Dict[str, str], Dict[str, str]]]):
    # Return directly if it's just one gene
    if len(genes) == 1:
        return genes[0][1]["name"]

    # For a complex, return str joined by " and "
    if genes[0][1]["label"] == "GeneProductComplex":
        return " and ".join([x[1]["name"] for x in genes[1:]])

    # GeneProductSet could contain genes or complexes
    res_complex, res_genes = {}, []
    for l, r in genes[1:]:
        # Everything in the set
        if l["label"] == "GeneProductSet":
            # GeneProductSet -> GeneProduct
            if r["label"] == "GeneProduct":
                res_genes.append(r["name"])
            # GeneProductSet -> GeneProductComplex
            else:
                res_complex[r["metaId"]] = []
        else:
            # GeneProductSet -> GeneProductComplex -> GeneProduct
            res_complex[l["metaId"]].append(r["name"])

    res_complex = {k: " and ".join(v) for k, v in res_complex.items()}
    rule_complex = [f"( {v} )" for v in res_complex.values()]

    res = " or ".join(res_genes + rule_complex)
    return res
