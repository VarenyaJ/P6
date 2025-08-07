"""
Measurement domain model.

Defines the MeasurementRecord dataclass for capturing quantitative measurements.
"""

from dataclasses import dataclass


@dataclass
class MeasurementRecord:
    """
    Represents a measurement entry for a patient.

    Attributes:
        patient_ID: Unique alphanumeric patient identifier.
        measurement_type: CURIE of measurement type (e.g. 'LOINC:4548-4').
        measurement_value: Numeric value of the measurement.
        measurement_unit: Unit CURIE or string (e.g. 'mmol/L').
        measurement_timestamp: ISO timestamp string (e.g. '2025-07-15T14:23:00').
    """

    patient_ID: str
    measurement_type: str
    measurement_value: float
    measurement_unit: str
    measurement_timestamp: str
