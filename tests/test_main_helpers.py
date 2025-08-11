"""
Unit tests for small helpers in __main__.py:
- _prepare_output_dir: creates timestamped folder
- _report_issues: prints warnings/errors to stdout
"""

import re
from P6.__main__ import _prepare_output_dir, _report_issues
from stairval.notepad import create_notepad


def test_prepare_output_dir_creates_timestamped_folder(tmp_path, monkeypatch):
    """
    The output path should look like:
    <repository rootH>/phenopacket_from_excel/YYYY-MM-DD_HH-MM-SS/phenopackets
    """
    monkeypatch.chdir(tmp_path)
    out = _prepare_output_dir()
    assert out.exists() and out.is_dir()
    assert re.search(r"\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}/phenopackets$", str(out))


def test_report_issues_outputs_both_blocks(capsys):
    """
    When notepad contains both warnings and errors, the helper should print both sections.
    """
    n = create_notepad("report")
    n.add_warning("warn 1")
    n.add_error("err 1")

    _report_issues(n)
    out = capsys.readouterr().out
    assert "Warnings found in mapping" in out
    assert "warn 1" in out
    assert "Errors found in mapping" in out
    assert "err 1" in out
