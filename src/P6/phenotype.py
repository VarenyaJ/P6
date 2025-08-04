"""
Phenotype domain model.

Defines the Phenotype class for HPO term annotations per patient visit.
"""

import re
from dataclasses import dataclass

# Patterns
_VALID_ID = re.compile(r"^[A-Za-z0-9]+$")
_HPO_ID_PATTERN = re.compile(r"^(?:HP:\d{7}|\d{7})$")
_TIMESTAMP_PATTERN = re.compile(r"^T\d+$")


@dataclass
class Phenotype:
    """
    Represents a single HPO annotation for a patient at a given timestamp.

    Attributes:
        phenotype_patient_ID: Unique alphanumeric patient identifier.
        HPO_ID: HPO term identifier (“HP:0001250” or “0001250”).
        date_of_observation: Visit timestamp (T0, T1, T2, …).
        status: True if observed, False if excluded.
    """

    phenotype_patient_ID: str
    HPO_ID: str
    date_of_observation: str
    status: bool

    def __post_init__(self):
        # Validate patient ID
        if not _VALID_ID.match(self.phenotype_patient_ID):
            raise ValueError(f"Invalid patient ID: {self.phenotype_patient_ID!r}")

        # Validate HPO ID
        if not _HPO_ID_PATTERN.match(self.HPO_ID):
            raise ValueError(f"Invalid HPO ID: {self.HPO_ID!r}")

        # Validate timestamp
        if not isinstance(
            self.date_of_observation, str
        ) or not _TIMESTAMP_PATTERN.match(self.date_of_observation):
            raise ValueError(
                f"Invalid date_of_observation: {self.date_of_observation!r}"
            )

        # Validate status
        if not isinstance(self.status, bool):
            raise ValueError(
                f"status must be a boolean, got {type(self.status).__name__}"
            )
