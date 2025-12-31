"""
HL7 v2 Message Parser

This module parses HL7 v2 messages (ORU^R01 - Lab Results) to extract
patient demographics and observation data for integration with RNA-Seq analysis.
"""

import logging
import re
from typing import Dict, List, Optional, Any
from datetime import datetime
from dataclasses import dataclass, field

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@dataclass
class PatientDemographics:
    """Patient demographic information."""
    patient_id: str
    patient_name: str
    date_of_birth: Optional[str] = None
    sex: Optional[str] = None
    address: Optional[str] = None
    phone: Optional[str] = None


@dataclass
class Observation:
    """HL7 observation (OBX segment)."""
    observation_id: str
    value_type: str
    observation_identifier: str
    observation_value: str
    units: Optional[str] = None
    reference_range: Optional[str] = None
    abnormal_flags: Optional[str] = None
    observation_status: str = "F"


@dataclass
class HL7Message:
    """Parsed HL7 v2 message."""
    message_type: str
    message_datetime: str
    patient: Optional[PatientDemographics] = None
    observations: List[Observation] = field(default_factory=list)
    raw_message: str = ""


class HL7v2Parser:
    """
    Parser for HL7 v2.x messages.
    
    Supports ORU^R01 (Observation Result/Unsolicited) messages commonly
    used for laboratory results.
    """
    
    # Field separators
    FIELD_SEP = '|'
    COMPONENT_SEP = '^'
    REPEAT_SEP = '~'
    ESCAPE_CHAR = '\\'
    SUBCOMPONENT_SEP = '&'
    
    def __init__(self):
        """Initialize HL7 parser."""
        logger.info("Initialized HL7 v2 parser")
    
    def parse_message(self, hl7_message: str) -> Optional[HL7Message]:
        """
        Parse HL7 v2 message into structured format.
        
        Args:
            hl7_message: Raw HL7 message string
            
        Returns:
            Parsed HL7Message object or None if error
        """
        try:
            logger.info("Parsing HL7 v2 message")
            
            # Split message into segments
            segments = [line.strip() for line in hl7_message.strip().split('\n') if line.strip()]
            
            if not segments:
                logger.error("Empty HL7 message")
                return None
            
            # Parse MSH (Message Header) segment
            msh_segment = segments[0]
            if not msh_segment.startswith('MSH'):
                logger.error("Invalid HL7 message: missing MSH segment")
                return None
            
            message = HL7Message(
                message_type="",
                message_datetime="",
                raw_message=hl7_message
            )
            
            # Parse each segment
            for segment in segments:
                segment_type = segment[:3]
                
                if segment_type == 'MSH':
                    self._parse_msh(segment, message)
                elif segment_type == 'PID':
                    self._parse_pid(segment, message)
                elif segment_type == 'OBX':
                    self._parse_obx(segment, message)
            
            logger.info(f"Parsed HL7 message: {message.message_type}")
            return message
            
        except Exception as e:
            logger.error(f"Error parsing HL7 message: {str(e)}")
            return None
    
    def _parse_msh(self, segment: str, message: HL7Message) -> None:
        """Parse MSH (Message Header) segment."""
        fields = segment.split(self.FIELD_SEP)
        
        if len(fields) < 9:
            logger.warning("Incomplete MSH segment")
            return
        
        # Message type (field 9)
        message_type_field = fields[8].split(self.COMPONENT_SEP)
        if message_type_field:
            message.message_type = message_type_field[0]
        
        # Message datetime (field 7)
        if len(fields) > 6:
            message.message_datetime = self._parse_datetime(fields[6])
    
    def _parse_pid(self, segment: str, message: HL7Message) -> None:
        """Parse PID (Patient Identification) segment."""
        fields = segment.split(self.FIELD_SEP)
        
        if len(fields) < 6:
            logger.warning("Incomplete PID segment")
            return
        
        # Patient ID (field 3)
        patient_id = fields[3] if len(fields) > 3 else ""
        if self.COMPONENT_SEP in patient_id:
            patient_id = patient_id.split(self.COMPONENT_SEP)[0]
        
        # Patient name (field 5)
        patient_name = ""
        if len(fields) > 5:
            name_components = fields[5].split(self.COMPONENT_SEP)
            if len(name_components) >= 2:
                # Last^First^Middle format
                patient_name = f"{name_components[1]} {name_components[0]}"
            else:
                patient_name = fields[5]
        
        # Date of birth (field 7)
        dob = self._parse_datetime(fields[7]) if len(fields) > 7 else None
        
        # Sex (field 8)
        sex = fields[8] if len(fields) > 8 else None
        
        # Address (field 11)
        address = None
        if len(fields) > 11:
            address_components = fields[11].split(self.COMPONENT_SEP)
            address = ', '.join([c for c in address_components if c])
        
        # Phone (field 13)
        phone = fields[13] if len(fields) > 13 else None
        
        message.patient = PatientDemographics(
            patient_id=patient_id,
            patient_name=patient_name,
            date_of_birth=dob,
            sex=sex,
            address=address,
            phone=phone
        )
    
    def _parse_obx(self, segment: str, message: HL7Message) -> None:
        """Parse OBX (Observation/Result) segment."""
        fields = segment.split(self.FIELD_SEP)
        
        if len(fields) < 6:
            logger.warning("Incomplete OBX segment")
            return
        
        # Set ID (field 1)
        observation_id = fields[1] if len(fields) > 1 else ""
        
        # Value type (field 2)
        value_type = fields[2] if len(fields) > 2 else ""
        
        # Observation identifier (field 3)
        obs_identifier_field = fields[3] if len(fields) > 3 else ""
        obs_identifier_parts = obs_identifier_field.split(self.COMPONENT_SEP)
        observation_identifier = obs_identifier_parts[0] if obs_identifier_parts else ""
        
        # Observation value (field 5)
        observation_value = fields[5] if len(fields) > 5 else ""
        
        # Units (field 6)
        units = None
        if len(fields) > 6:
            units_field = fields[6].split(self.COMPONENT_SEP)
            units = units_field[0] if units_field else None
        
        # Reference range (field 7)
        reference_range = fields[7] if len(fields) > 7 else None
        
        # Abnormal flags (field 8)
        abnormal_flags = fields[8] if len(fields) > 8 else None
        
        # Observation status (field 11)
        status = fields[11] if len(fields) > 11 else "F"
        
        observation = Observation(
            observation_id=observation_id,
            value_type=value_type,
            observation_identifier=observation_identifier,
            observation_value=observation_value,
            units=units,
            reference_range=reference_range,
            abnormal_flags=abnormal_flags,
            observation_status=status
        )
        
        message.observations.append(observation)
    
    def _parse_datetime(self, dt_string: str) -> str:
        """
        Parse HL7 datetime format to ISO format.
        
        Args:
            dt_string: HL7 datetime string (e.g., "20240101120000")
            
        Returns:
            ISO format datetime string
        """
        try:
            # HL7 format: YYYYMMDDHHMMSS
            if len(dt_string) >= 8:
                year = dt_string[0:4]
                month = dt_string[4:6]
                day = dt_string[6:8]
                
                hour = dt_string[8:10] if len(dt_string) >= 10 else "00"
                minute = dt_string[10:12] if len(dt_string) >= 12 else "00"
                second = dt_string[12:14] if len(dt_string) >= 14 else "00"
                
                dt = datetime(
                    int(year), int(month), int(day),
                    int(hour), int(minute), int(second)
                )
                return dt.isoformat()
            
            return dt_string
            
        except Exception as e:
            logger.warning(f"Error parsing datetime: {str(e)}")
            return dt_string
    
    def to_dict(self, message: HL7Message) -> Dict[str, Any]:
        """
        Convert parsed HL7 message to dictionary.
        
        Args:
            message: Parsed HL7Message
            
        Returns:
            Dictionary representation
        """
        result = {
            'message_type': message.message_type,
            'message_datetime': message.message_datetime,
            'patient': None,
            'observations': []
        }
        
        if message.patient:
            result['patient'] = {
                'patient_id': message.patient.patient_id,
                'patient_name': message.patient.patient_name,
                'date_of_birth': message.patient.date_of_birth,
                'sex': message.patient.sex,
                'address': message.patient.address,
                'phone': message.patient.phone
            }
        
        for obs in message.observations:
            result['observations'].append({
                'observation_id': obs.observation_id,
                'value_type': obs.value_type,
                'observation_identifier': obs.observation_identifier,
                'observation_value': obs.observation_value,
                'units': obs.units,
                'reference_range': obs.reference_range,
                'abnormal_flags': obs.abnormal_flags,
                'observation_status': obs.observation_status
            })
        
        return result
    
    def map_to_rnaseq_metadata(self, message: HL7Message) -> Dict[str, Any]:
        """
        Map HL7 message to RNA-Seq metadata format.
        
        Args:
            message: Parsed HL7Message
            
        Returns:
            RNA-Seq metadata dictionary
        """
        metadata = {
            'sample_id': message.patient.patient_id if message.patient else '',
            'collection_date': message.message_datetime,
            'patient_demographics': {}
        }
        
        if message.patient:
            metadata['patient_demographics'] = {
                'patient_id': message.patient.patient_id,
                'age': self._calculate_age(message.patient.date_of_birth),
                'sex': message.patient.sex
            }
        
        # Map relevant observations to metadata
        metadata['clinical_observations'] = []
        for obs in message.observations:
            metadata['clinical_observations'].append({
                'test_code': obs.observation_identifier,
                'test_value': obs.observation_value,
                'units': obs.units,
                'reference_range': obs.reference_range,
                'abnormal': obs.abnormal_flags
            })
        
        return metadata
    
    def _calculate_age(self, dob: Optional[str]) -> Optional[int]:
        """Calculate age from date of birth."""
        if not dob:
            return None
        
        try:
            dob_dt = datetime.fromisoformat(dob)
            today = datetime.now()
            age = today.year - dob_dt.year
            if today.month < dob_dt.month or (today.month == dob_dt.month and today.day < dob_dt.day):
                age -= 1
            return age
        except:
            return None


# Example usage
if __name__ == "__main__":
    # Sample HL7 v2 ORU^R01 message
    sample_hl7_message = """MSH|^~\\&|LAB|HOSPITAL|RIS|HOSPITAL|20240101120000||ORU^R01|MSG123456|P|2.5
PID|1||PAT123456^^^HOSPITAL^MR||DOE^JOHN^A||19800115|M|||123 MAIN ST^^CITY^STATE^12345^USA|(555)123-4567
OBR|1|ORD123456|RES123456|PANEL^RNASEQ PANEL^L|||20240101100000
OBX|1|NM|GENE001^TP53 Expression^L||1234.56|TPM|100-2000||||F
OBX|2|NM|GENE002^BRCA1 Expression^L||567.89|TPM|50-1000||||F
OBX|3|NM|GENE003^EGFR Expression^L||2345.67|TPM|200-3000|H|||F"""
    
    # Parse message
    parser = HL7v2Parser()
    parsed_message = parser.parse_message(sample_hl7_message)
    
    if parsed_message:
        print("=== Parsed HL7 Message ===")
        import json
        print(json.dumps(parser.to_dict(parsed_message), indent=2))
        print()
        
        print("=== RNA-Seq Metadata ===")
        metadata = parser.map_to_rnaseq_metadata(parsed_message)
        print(json.dumps(metadata, indent=2))
    else:
        print("Failed to parse HL7 message")
