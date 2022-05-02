import libsbml
import pytest
from libsbml import Model, SBMLDocument
from metabolike.parser import SBMLParser


@pytest.fixture
def sbml_parser():
    # No need to pass an actual database to the parser
    with pytest.raises(ValueError):
        SBMLParser(None)
    return SBMLParser("tests/data/metabolic-reactions.sbml")


@pytest.fixture
def sbml_model(sbml_parser: SBMLParser):
    doc = sbml_parser.read_sbml(sbml_parser.sbml_file)
    assert isinstance(doc, SBMLDocument)

    model = doc.getModel()
    assert isinstance(model, Model)

    return model


def test_collect_compartments(sbml_parser: SBMLParser, sbml_model: Model):
    assert sbml_model.getNumCompartments() == 33
    compartments = sbml_model.getListOfCompartments()
    nodes = sbml_parser.collect_compartments(compartments)

    assert isinstance(nodes, list)
    assert len(nodes) == 33
    assert all("metaId" in n for n in nodes)
    assert all("props" in n for n in nodes)
    assert nodes[-1] == {"metaId": "x", "props": {"name": "CCO-PEROX-LUM"}}


def test_collect_compounds(sbml_parser: SBMLParser, sbml_model: Model):
    assert sbml_model.getNumSpecies() == 45
    species = sbml_model.getListOfSpecies()
    nodes = sbml_parser.collect_compounds(species)

    assert isinstance(nodes, list)
    assert len(nodes) == 45
    assert all("metaId" in n for n in nodes)
    assert all("props" in n for n in nodes)
    assert all("compartment" in n for n in nodes)
    assert nodes[0] == {
        "metaId": "G3P_c",
        "props": {
            "name": "3-phospho-D-glycerate",
            "chemicalFormula": "C3H4O7P",
            "boundaryCondition": False,
            "hasOnlySubstanceUnits": False,
            "constant": False,
        },
        "compartment": "c",
        "rdf": [
            {
                "bioQual": "is",
                "rdf": {
                    "biocyc": "G3P",
                    "inchikey": "OSJPPGNTCRNQQC-UWTATZPHSA-K",
                    "chebi": "58272",
                    "cas": "820-11-1",
                    "hmdb": "HMDB60180",
                    "pubchemCompound": "25245548",
                    "keggCompound": "C00197",
                    "metabolights": "MTBLC58272",
                },
            },
            {"bioQual": "hasProperty", "rdf": {"sbo": "0000247"}},
        ],
    }


def test_collect_gene_products(sbml_parser: SBMLParser, sbml_model: Model):
    fbc = sbml_model.getPlugin("fbc")
    assert fbc.getNumGeneProducts() == 118
    gene_products = fbc.getListOfGeneProducts()
    nodes = sbml_parser.collect_gene_products(gene_products)

    assert isinstance(nodes, list)
    assert len(nodes) == 118
    assert all("metaId" in n for n in nodes)
    assert all("props" in n for n in nodes)
    assert all("rdf" in n for n in nodes)

    assert nodes[0] == {
        "metaId": "G_G_3784",
        "props": {"name": "apgM", "label": "G-3784"},
        "rdf": [
            {
                "bioQual": "is",
                "rdf": {
                    "uniprot": "Q9X295",
                },
            },
            {
                "bioQual": "isEncodedBy",
                "rdf": {"biocyc": "G-3784", "ncbigene": "897849"},
            },
        ],
    }


def test_collect_reactions(sbml_parser: SBMLParser, sbml_model: Model):
    assert sbml_model.getNumReactions() == 12
    reactions = sbml_model.getListOfReactions()
    nodes = sbml_parser.collect_reactions(reactions)

    assert isinstance(nodes, list)
    assert len(nodes) == 12
    assert all("metaId" in n for n in nodes)
    assert all("props" in n for n in nodes)
    assert all("rdf" in n for n in nodes)
    assert all("reactants" in n for n in nodes)
    assert all("products" in n for n in nodes)

    assert nodes[0] == {
        "metaId": "RXN__45__15513",
        "props": {"name": "RXN-15513", "fast": False, "reversible": True},
        "rdf": [
            {"bioQual": "is", "rdf": {"biocyc": "RXN-15513", "ecCode": ["5.4.2.11"]}}
        ],
        "reactants": [
            {"cpdId": "_2__45__PG_c", "props": {"stoichiometry": 1.0, "constant": True}}
        ],
        "products": [
            {"cpdId": "G3P_c", "props": {"stoichiometry": 1.0, "constant": True}}
        ],
    }


@pytest.fixture
def complex_node():
    n = libsbml.FbcAnd()
    g1: libsbml.GeneProductRef = n.createGeneProductRef()
    g1.setGeneProduct("g1")

    g2: libsbml.FbcOr = n.createOr()
    g2_1: libsbml.GeneProductRef = g2.createGeneProductRef()
    g2_1.setGeneProduct("g2_1")
    g2_2: libsbml.GeneProductRef = g2.createGeneProductRef()
    g2_2.setGeneProduct("g2_2")
    return n


@pytest.fixture
def set_node():
    n = libsbml.FbcOr()
    g1: libsbml.GeneProductRef = n.createGeneProductRef()
    g1.setGeneProduct("gp1")

    g2: libsbml.FbcAnd = n.createAnd()
    g2_1: libsbml.GeneProductRef = g2.createGeneProductRef()
    g2_1.setGeneProduct("gp2_1")
    g2_2: libsbml.GeneProductRef = g2.createGeneProductRef()
    g2_2.setGeneProduct("gp2_2")
    return n


def test_collect_reaction_gene_product_links(
    sbml_parser: SBMLParser, sbml_model: Model
):
    reactions = sbml_model.getListOfReactions()
    (
        rxn_genes,
        gene_sets,
        gene_complexes,
    ) = sbml_parser.collect_reaction_gene_product_links(reactions)

    assert isinstance(rxn_genes, dict)
    num_rxns = len(rxn_genes)
    for i, v in zip(range(num_rxns), rxn_genes.values()):
        assert v == f"GeneProductSet{i}"

    assert len(gene_sets) == 12
    assert len(gene_sets["GeneProductSet0"]) == 6

    assert len(gene_complexes) == 1
    assert gene_complexes["GeneProductComplex0"] == {"G_YGR240C", "G_YMR205C"}


def test_parse_gene_association_nodes(
    sbml_parser: SBMLParser, complex_node: libsbml.FbcAnd, set_node: libsbml.FbcOr
):
    gene_sets, gene_complexes = {}, {}
    n = libsbml.GeneProductRef()
    n.setGeneProduct("gene1")
    res = sbml_parser._parse_gene_association_nodes(n, gene_sets, gene_complexes)
    assert res == "gene1"
    assert gene_sets == {}
    assert gene_complexes == {}

    res = sbml_parser._parse_gene_association_nodes(
        complex_node, gene_sets, gene_complexes
    )
    assert res == "GeneProductComplex0"
    assert len(gene_sets) == 1
    assert len(gene_complexes) == 1
    assert res in gene_complexes

    res = sbml_parser._parse_gene_association_nodes(set_node, gene_sets, gene_complexes)
    assert res == "GeneProductSet1"
    assert len(gene_sets) == 2
    assert len(gene_complexes) == 2
    assert res in gene_sets

    with pytest.raises(ValueError):
        sbml_parser._parse_gene_association_nodes(
            libsbml.Reaction(3, 1), gene_sets, gene_complexes
        )


def test_parse_gene_product_complex(
    sbml_parser: SBMLParser, complex_node: libsbml.FbcAnd
):
    gene_sets, gene_complexes = {}, {}

    assert (
        sbml_parser._parse_gene_product_complex(
            libsbml.FbcAnd(), gene_sets, gene_complexes
        )
        is None
    )

    k = sbml_parser._parse_gene_product_complex(complex_node, gene_sets, gene_complexes)
    assert k == "GeneProductComplex0"

    assert len(gene_complexes) == 1
    assert gene_complexes[k] == {"g1", "GeneProductSet0"}

    assert len(gene_sets) == 1
    assert gene_sets["GeneProductSet0"] == {"g2_1", "g2_2"}

    with pytest.raises(ValueError):
        complex_node.createAnd()
        sbml_parser._parse_gene_product_complex(complex_node, gene_sets, gene_complexes)


def test_parse_gene_product_set(sbml_parser: SBMLParser, set_node: libsbml.FbcOr):
    gene_sets, gene_complexes = {}, {}
    assert (
        sbml_parser._parse_gene_product_set(libsbml.FbcOr(), gene_sets, gene_complexes)
        is None
    )

    k = sbml_parser._parse_gene_product_set(set_node, gene_sets, gene_complexes)
    assert k == "GeneProductSet0"

    assert len(gene_sets) == 1
    assert gene_sets[k] == {"gp1", "GeneProductComplex0"}

    assert len(gene_complexes) == 1
    assert gene_complexes["GeneProductComplex0"] == {"gp2_1", "gp2_2"}

    with pytest.raises(ValueError):
        set_node.createOr()
        sbml_parser._parse_gene_product_set(set_node, gene_sets, gene_complexes)


def test_add_gene_product_group(sbml_parser: SBMLParser):
    g = {}
    k0 = sbml_parser._add_gene_product_group({"g1", "g2"}, g, "GeneSet")
    assert k0 == "GeneSet0"
    assert g == {k0: {"g1", "g2"}}

    k1 = sbml_parser._add_gene_product_group({"g2", "g3"}, g, "GeneSet")
    assert len(g) == 2
    assert k1 == "GeneSet1"
    assert g[k1] == {"g2", "g3"}

    k2 = sbml_parser._add_gene_product_group({"g2", "g1"}, g, "GeneSet")
    assert len(g) == 2
    assert k2 == k0


def test_split_uri(sbml_parser: SBMLParser):
    assert sbml_parser._split_uri("http://identifiers.org/ec-code/1.1.1.1") == (
        "ec-code",
        "1.1.1.1",
    )
    assert sbml_parser._split_uri("http://identifiers.org/biocyc/META:RXN-15513") == (
        "biocyc",
        "RXN-15513",
    )
    assert sbml_parser._split_uri("http://identifiers.org/kegg.reaction/R00199") == (
        "kegg-reaction",
        "R00199",
    )
