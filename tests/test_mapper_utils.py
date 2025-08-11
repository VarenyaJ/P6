"""
Small utility tests for mapper statics:
- _normalize_time_like
- _to_bool
"""

from P6.mapper import DefaultMapper


def test_normalize_time_like_variants():
    n = DefaultMapper._normalize_time_like
    assert n(20200101) == "T20200101"
    assert n("T2020") == "T2020"
    assert n(" 2020 ") == "T2020"
    assert n("") == ""
    assert n(None) == ""


def test_to_bool_truth_table():
    b = DefaultMapper._to_bool
    for t in [1, "1", "true", "TRUE", "Yes", "y", True]:
        assert b(t) is True
    for f in [0, "0", "false", "no", "", None, False]:
        assert b(f) is False
