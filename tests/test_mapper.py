import pytest

import hpotk
import pandas as pd

from stairval.notepad import create_notepad

from P6.mapper import DefaultMapper


class TestDefaultMapper:
    @pytest.fixture(scope="class")
    def mapper(self, hpo: hpotk.MinimalOntology) -> DefaultMapper:
        return DefaultMapper(hpo)

    def test_apply_mapping(self, mapper: DefaultMapper):
        notepad = create_notepad("whatever")

        pps = mapper.apply_mapping(
            {},  # TODO:
            notepad,
        )

        assert len(pps) is not None

