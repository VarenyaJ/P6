import gzip
import hpotk
import os
import pytest
import shutil

from pathlib import Path


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

# use the `hpo` already defined above in this file
@pytest.fixture
def decompressed_hpo(fpath_hpo: str, tmp_path: Path) -> str:
    """
    Decompress the gzipped HPO JSON to a plain JSON file and return its path.
    """
    out = tmp_path / "hp.json"
    with gzip.open(fpath_hpo, "rb") as f_in, open(out, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    return str(out)