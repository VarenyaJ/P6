import os
import pytest

import hpotk


@pytest.fixture(scope="session")
def fpath_test_dir() -> str:
    """
    Path to `tests` folder.
    """
    return os.path.dirname(os.path.abspath(__file__))


@pytest.fixture(scope="session")
def fpath_hpo(fpath_test_dir: str) -> str:
    return os.path.join(fpath_test_dir, "data", "hp.v2024-04-26.json.gz")


@pytest.fixture(scope="session")
def hpo(fpath_hpo: str) -> hpotk.MinimalOntology:
    return hpotk.load_minimal_ontology(fpath_hpo)
