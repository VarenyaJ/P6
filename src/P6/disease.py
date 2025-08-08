"""
Disease domain model.

Defines the DiseaseRecord dataclass for capturing disease annotations.
"""

from dataclasses import dataclass


@dataclass
class DiseaseRecord:
    """
    Represents a disease entry for a patient.

    Attributes:
        patient_ID: Unique alphanumeric patient identifier.
        disease_term: CURIE of the disease term (e.g. 'OMIM:266600').
        disease_label: Human-readable label for the disease.
        disease_onset: Date string in 'YYYY-MM-DD' format.
        disease_status: True if the disease is present, False if excluded.
    """

    patient_ID: str
    disease_term: str
    disease_label: str
    disease_onset: str
    disease_status: bool
