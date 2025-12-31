"""
Test healthcare standards modules
"""
import pytest
from healthcare_standards.fhir_genomics_integration import create_genomic_observation
from healthcare_standards.hl7_v2_parser import parse_hl7_message, create_hl7_message
from healthcare_standards.omop_cdm_mapping import map_to_omop


def test_fhir_genomics():
    """Test FHIR genomics integration"""
    result = create_genomic_observation("patient-123", {"BRCA1": 234.5, "TP53": 156.2})
    assert result is not None
    assert result["resourceType"] == "Observation"
    assert len(result["component"]) == 2


def test_hl7_parser():
    """Test HL7 v2 parser"""
    msg = "MSH|^~\\&|LAB|HOSP|||20231231120000||ORU^R01|001|P|2.5"
    result = parse_hl7_message(msg)
    assert result is not None
    assert result["message_type"] == "ORU^R01"


def test_hl7_creator():
    """Test HL7 message creation"""
    msg = create_hl7_message("patient-123", [
        {"test": "GLU", "value": 95, "unit": "mg/dL"}
    ])
    assert msg is not None
    assert "MSH" in msg
    assert "patient-123" in msg


def test_omop_mapping():
    """Test OMOP CDM mapping"""
    data = {
        "patient_id": "123",
        "lab_values": {"glucose": 95.0},
        "genomic_data": {"BRCA1": 234.5}
    }
    result = map_to_omop(data)
    assert result is not None
    assert "person" in result
    assert len(result["measurement"]) == 1
    assert len(result["observation"]) == 1
