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
