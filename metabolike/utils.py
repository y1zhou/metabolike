from typing import Dict, List


def generate_gene_reaction_rule(genes: List[List[Dict[str, str]]]):
    # Return directly if it's just one gene
    if len(genes) == 1:
        return genes[0][1]["name"]

    # For a complex, return str joined by " and "
    if genes[0][1]["label"] == "Complex":
        return " and ".join([x[1]["name"] for x in genes[1:]])

    # EntitySet could contain genes or complexes
    res_complex, res_genes = {}, []
    for l, r in genes[1:]:
        # Everything in the set
        if l["label"] == "EntitySet":
            # EntitySet -> GeneProduct
            if r["label"] == "GeneProduct":
                res_genes.append(r["name"])
            # EntitySet -> Complex
            else:
                res_complex[r["mcId"]] = []
        else:
            # EntitySet -> Complex -> GeneProduct
            res_complex[l["mcId"]].append(r["name"])

    res_complex = {k: " and ".join(v) for k, v in res_complex.items()}
    rule_complex = [f"( {v} )" for v in res_complex.values()]

    res = " or ".join(res_genes + rule_complex)
    return res
