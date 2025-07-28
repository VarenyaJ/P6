import os
import re
import pytest
from click.testing import CliRunner
from P6.scripts.__main__ import main

@pytest.fixture
def sample_xlsx_path():
    # adjust this if your test files live elsewhere
    return os.path.join(os.path.dirname(__file__), "..", "Python_headers_phenocopy_transformation.xlsx")

def test_parse_excel_creates_records(sample_xlsx_path):
    """
    Runs `p6 parse-excel` on the sample workbook and
    asserts that it reports nonzero Genotype & Phenotype counts.
    """
    runner = CliRunner()
    result = runner.invoke(main, ["parse-excel", sample_xlsx_path])

    # it should exit cleanly
    assert result.exit_code == 0, result.output

    # output should contain "Created N Genotype objects"
    m1 = re.search(r"Created (\d+) Genotype objects", result.output)
    assert m1, f"Missing Genotype line in output:\n{result.output}"
    assert int(m1.group(1)) > 0, "Expected at least one Genotype record"

    # and similarly for Phenotype
    m2 = re.search(r"Created (\d+) Phenotype objects", result.output)
    assert m2, f"Missing Phenotype line in output:\n{result.output}"
    assert int(m2.group(1)) > 0, "Expected at least one Phenotype record"
