import glob
import gzip
import os
import pytest
import re
import tempfile
from click.testing import CliRunner
from P6.__main__ import main

import pytest
import gzip
import tempfile
import os

@pytest.fixture(scope="session", autouse=True)
def verify_hpo_file():
    # Ensure HPO file exists before running any tests
    data_dir = os.path.join(os.path.dirname(__file__), "data")
    hpo_path = os.path.join(data_dir, 'hp.v2024-04-26.json.gz')
    assert os.path.exists(hpo_path), f"HPO file not found: {hpo_path}"


@pytest.fixture(scope="function")
def decompressed_hpo():
    # Provides temporary decompressed HPO file for testing
    data_dir = os.path.join(os.path.dirname(__file__), "data")
    compressed_path = os.path.join(data_dir, 'hp.v2024-04-26.json.gz')

    with tempfile.NamedTemporaryFile(mode='w', delete=True) as temp_file:
        # Decompress to temporary file
        with gzip.open(compressed_path, 'rt') as gz_file:
            temp_file.write(gz_file.read())
        temp_file.flush()

        yield temp_file.name

        # File should get deleted automatically

@pytest.mark.parametrize("sample_xlsx_path", glob.glob(os.path.join(os.path.dirname(__file__), "data", "*.xlsx")),)
def test_parse_excel_creates_records(sample_xlsx_path, decompressed_hpo):
    """
    Runs `p6 parse-excel` on each workbook in tests/data/ and
    asserts that it reports nonzero Genotype & Phenotype counts.
    """
    runner = CliRunner()
    data_dir = os.path.join(os.path.dirname(__file__), "data")

    # Use the decompressed HPO file path in the command
    result = runner.invoke(main, ["parse-excel", sample_xlsx_path, "--d", data_dir, "--hpo", decompressed_hpo])
    # it should exit cleanly
    assert result.exit_code == 0, result.output

    # Print debug information
    print("\nDebug Output:")
    print(f"Full output:\n{result.output}")
    # Check for warnings about ambiguous sheets
    if "Skipping sheet: cannot unambiguously classify" in result.output:
        print("\nWarning: Some sheets were skipped due to ambiguous classification!")

    # output should contain "Created N Genotype objects"
    m1 = re.search(r"Created (\d+) Genotype objects", result.output)
    assert m1, f"Missing Genotype line in output:\n{result.output}"
    assert int(m1.group(1)) > 0, "Expected at least one Genotype record"

    # and similarly for Phenotype
    m2 = re.search(r"Created (\d+) Phenotype objects", result.output)
    assert m2, f"Missing Phenotype line in output:\n{result.output}"
    assert int(m2.group(1)) > 0, "Expected at least one Phenotype record"
