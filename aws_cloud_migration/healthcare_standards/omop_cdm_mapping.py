"""
OMOP Common Data Model (CDM) Mapping Module
Map clinical and genomic data to OMOP CDM standard
"""


def map_to_omop(data):
    """
    Map data to OMOP CDM format
    
    Args:
        data: Dictionary containing clinical/genomic data
        
    Returns:
        Dictionary in OMOP CDM format
    """
    omop_data = {
        "person": {
            "person_id": data.get("patient_id"),
            "gender_concept_id": 8507,  # Male (example)
            "year_of_birth": 1980,
            "race_concept_id": 8527,  # White (example)
            "ethnicity_concept_id": 38003564  # Not Hispanic (example)
        },
        "measurement": [],
        "observation": []
    }
    
    # Map lab values to measurements
    if "lab_values" in data:
        for lab_name, lab_value in data["lab_values"].items():
            measurement = {
                "measurement_id": None,  # Auto-generated
                "person_id": data.get("patient_id"),
                "measurement_concept_id": None,  # Map lab name to concept
                "measurement_date": None,
                "value_as_number": lab_value,
                "unit_concept_id": None
            }
            omop_data["measurement"].append(measurement)
    
    # Map genomic data to observations
    if "genomic_data" in data:
        for gene, expression in data["genomic_data"].items():
            observation = {
                "observation_id": None,
                "person_id": data.get("patient_id"),
                "observation_concept_id": None,  # Map gene to concept
                "observation_date": None,
                "value_as_number": expression
            }
            omop_data["observation"].append(observation)
    
    return omop_data


if __name__ == "__main__":
    # Example usage
    test_data = {
        "patient_id": "123",
        "lab_values": {"glucose": 95.0, "hemoglobin": 14.5},
        "genomic_data": {"BRCA1": 234.5}
    }
    result = map_to_omop(test_data)
    print(f"Mapped to OMOP with {len(result['measurement'])} measurements")
