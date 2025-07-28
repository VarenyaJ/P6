import pytest
from P6.periodicity import FrequencyModifier, Periodicity


def test_frequency_modifier_from_label():
    """Ensure that known labels map to the correct enum."""
    assert FrequencyModifier.from_label("Obligate") == FrequencyModifier.OBLIGATE
    assert (
        FrequencyModifier.from_label("very frequent") == FrequencyModifier.VERY_FREQUENT
    )


def test_frequency_modifier_invalid_label_raises():
    """Unknown labels must trigger a ValueError."""
    with pytest.raises(ValueError):
        FrequencyModifier.from_label("Sometimes")


def test_periodicity_wrapper_stores_enum():
    """Periodicity should accept a FrequencyModifier."""
    periodicity = Periodicity(frequency_modifier=FrequencyModifier.FREQUENT)
    assert isinstance(periodicity.frequency_modifier, FrequencyModifier)
