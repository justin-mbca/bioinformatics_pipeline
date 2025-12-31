"""
HL7 v2 Parser Module
Parse and create HL7 v2 messages for laboratory results
"""
from datetime import datetime


def parse_hl7_message(message):
    """
    Parse an HL7 v2 message
    
    Args:
        message: HL7 v2 message string
        
    Returns:
        Dictionary with parsed message components
    """
    lines = message.split('\n')
    parsed = {
        "message_type": None,
        "segments": []
    }
    
    for line in lines:
        if line.startswith("MSH"):
            parts = line.split('|')
            if len(parts) >= 9:
                parsed["message_type"] = parts[8]
        
        segment_type = line[:3] if len(line) >= 3 else ""
        parsed["segments"].append({
            "type": segment_type,
            "content": line
        })
    
    return parsed


def create_hl7_message(patient_id, lab_results):
    """
    Create an HL7 v2 ORU^R01 message for laboratory results
    
    Args:
        patient_id: Patient identifier
        lab_results: List of lab result dictionaries
        
    Returns:
        HL7 v2 message string
    """
    timestamp = datetime.now().strftime("%Y%m%d%H%M%S")
    message_id = f"MSG{timestamp}"
    
    # MSH segment
    msh = f"MSH|^~\\&|LAB|HOSPITAL|||{timestamp}||ORU^R01|{message_id}|P|2.5"
    
    # PID segment
    pid = f"PID|1||{patient_id}||DOE^JOHN||19800101|M"
    
    # OBR segment
    obr = f"OBR|1|{message_id}|{message_id}|LAB^Laboratory Panel|||{timestamp}"
    
    # OBX segments for each result
    obx_segments = []
    for i, result in enumerate(lab_results, 1):
        test_name = result.get("test", "UNKNOWN")
        value = result.get("value", "")
        unit = result.get("unit", "")
        obx = f"OBX|{i}|NM|{test_name}^{test_name}||{value}|{unit}|||N"
        obx_segments.append(obx)
    
    message = "\n".join([msh, pid, obr] + obx_segments)
    return message


if __name__ == "__main__":
    # Example usage
    msg = "MSH|^~\\&|LAB|HOSP|||20231231120000||ORU^R01|001|P|2.5"
    result = parse_hl7_message(msg)
    print(f"Parsed message type: {result['message_type']}")
    
    new_msg = create_hl7_message("patient-123", [
        {"test": "GLU", "value": 95, "unit": "mg/dL"},
        {"test": "HGB", "value": 14.5, "unit": "g/dL"}
    ])
    print(f"Created HL7 message: {len(new_msg.split(chr(10)))} segments")
