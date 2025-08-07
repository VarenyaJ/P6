"""
Biosample domain model.

Defines the BiosampleRecord dataclass for capturing sample metadata.
"""

from dataclasses import dataclass


@dataclass
class BiosampleRecord:
    """
    Represents a biosample entry for a patient.

    Attributes:
        patient_ID: Unique alphanumeric patient identifier.
        biosample_id: Unique identifier for the biosample.
        biosample_type: CURIE of the tissue or sample type (e.g. 'UBERON:0002107').
        collection_date: Date string in 'YYYY-MM-DD' format.
    """

    patient_ID: str
    biosample_id: str
    biosample_type: str
    collection_date: str
