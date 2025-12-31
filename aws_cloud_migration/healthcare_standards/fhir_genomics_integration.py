"""
FHIR Genomics Integration Module
Converts genomic data to FHIR format
"""

def create_genomic_observation(patient_id, gene_expressions):
    """
    Create a FHIR Observation resource for genomic data
    
    Args:
        patient_id: Patient identifier
        gene_expressions: Dictionary of gene name to expression value
        
    Returns:
        FHIR Observation resource as dictionary
    """
    observation = {
        "resourceType": "Observation",
        "id": f"genomic-{patient_id}",
        "status": "final",
        "category": [{
            "coding": [{
                "system": "http://terminology.hl7.org/CodeSystem/observation-category",
                "code": "laboratory",
                "display": "Laboratory"
            }]
        }],
        "code": {
            "coding": [{
                "system": "http://loinc.org",
                "code": "81247-9",
                "display": "Master HL7 genetic variant reporting panel"
            }]
        },
        "subject": {
            "reference": f"Patient/{patient_id}"
        },
        "component": []
    }
    
    for gene_name, expression_value in gene_expressions.items():
        component = {
            "code": {
                "coding": [{
                    "system": "http://www.genenames.org",
                    "code": gene_name,
                    "display": gene_name
                }]
            },
            "valueQuantity": {
                "value": expression_value,
                "unit": "TPM",
                "system": "http://unitsofmeasure.org",
                "code": "1"
            }
        }
        observation["component"].append(component)
    
    return observation


if __name__ == "__main__":
    # Example usage
    result = create_genomic_observation("patient-123", {"BRCA1": 234.5, "TP53": 156.2})
    print(f"Created FHIR observation with {len(result['component'])} gene components")
