"""
Periodicity domain model.

Defines frequency modifiers for cohort‑level phenotypic occurrence rates.
"""

from dataclasses import dataclass
from enum import Enum, auto


class FrequencyModifier(Enum):
    """
    Enumeration of frequency modifiers for phenotypic abnormalities.
    Mirrors HP_0040280–HP_0040285 under the “Frequency” class in HPO.
    """
    OBLIGATE = auto()
    VERY_FREQUENT = auto()
    FREQUENT = auto()
    OCCASIONAL = auto()
    VERY_RARE = auto()
    EXCLUDED = auto()

    @classmethod
    def from_label(cls, label: str) -> "FrequencyModifier":
        """
        Convert a human‑readable label into the corresponding enum.
        Strips punctuation and normalizes spacing/casing.
        """
        key = label.strip().lower().replace(" ", "_").replace("(", "").replace(")", "")
        mapping = {
            "obligate": cls.OBLIGATE,
            "very_frequent": cls.VERY_FREQUENT,
            "frequent": cls.FREQUENT,
            "occasional": cls.OCCASIONAL,
            "very_rare": cls.VERY_RARE,
            "excluded": cls.EXCLUDED,
        }
        try:
            return mapping[key]
        except KeyError:
            raise ValueError(f"Unknown frequency modifier label: {label!r}")


@dataclass
class Periodicity:
    """
    Wraps a FrequencyModifier for convenience in data models.
    """
    frequency_modifier: FrequencyModifier
