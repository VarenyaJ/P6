import hpotk
import os
import pytest


@pytest.fixture(scope="session")
def fpath_test_dir() -> str:
    """
    Path to `tests/data/` folder.
    """
    return os.path.join(os.path.dirname(__file__), "data")


@pytest.fixture(scope="session")
def fpath_hpo(fpath_test_dir: str) -> str:
    return os.path.join(fpath_test_dir, "hp.v2024-04-26.json.gz")


@pytest.fixture(scope="session")
def hpo(fpath_hpo: str) -> hpotk.MinimalOntology:
    """
    The PATH to a JSON file of HPO terms and IDs.
    `hpotk` should be able to read this directly without manual decompression.
    """
    return hpotk.load_minimal_ontology(fpath_hpo)
