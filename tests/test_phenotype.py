import pytest
from P6.phenotype import Phenotype


def test_valid_phenotype_instantiation():
    """A valid Phenotype should be created without error."""
    p = Phenotype(
        phenotype_patient_ID="PXYZ789",
        HPO_ID="HP:0001250",
        date_of_observation="T0",
        status=True
    )
    assert isinstance(p, Phenotype)


@pytest.mark.parametrize("bad_hpo", ["HP:123", "000123", "HP:ABCDEF1"])
def test_invalid_hpo_id_raises(bad_hpo):
    """Malformed HPO IDs must trigger a ValueError."""
    with pytest.raises(ValueError):
        Phenotype(
            phenotype_patient_ID="P1",
            HPO_ID=bad_hpo,
            date_of_observation="T1",
            status=False
        )


@pytest.mark.parametrize("bad_timestamp", ["0", "Time1", "T"])
def test_invalid_timestamp_raises(bad_timestamp):
    """Timestamps must be in the form T<integer>."""
    with pytest.raises(ValueError):
        Phenotype(
            phenotype_patient_ID="P2",
            HPO_ID="0001250",
            date_of_observation=bad_timestamp,
            status=True
        )