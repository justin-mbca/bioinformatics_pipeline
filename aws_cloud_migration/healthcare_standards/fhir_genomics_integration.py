"""
FHIR R4 Genomics Integration

This module provides integration with FHIR R4 (Fast Healthcare Interoperability Resources)
for representing genomic data, specifically RNA-Seq gene expression values.
"""

import json
import logging
from typing import Dict, List, Optional, Any
from datetime import datetime
from dataclasses import dataclass, asdict

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@dataclass
class FHIRCoding:
    """FHIR Coding data class."""
    system: str
    code: str
    display: str


@dataclass
class FHIRQuantity:
    """FHIR Quantity data class."""
    value: float
    unit: str
    system: str
    code: str


class FHIRGenomicsObservation:
    """
    Create FHIR R4 Observation resources for genomic data.
    
    Maps gene expression values to FHIR Observation resources with proper
    LOINC codes and HGNC gene identifiers.
    """
    
    # LOINC codes for genomic observations
    LOINC_GENE_EXPRESSION = "81247-9"  # Gene expression
    LOINC_GENE_ID = "48018-6"  # Gene studied [ID]
    LOINC_TPM = "LA6683-2"  # Transcripts per million (TPM)
    LOINC_FPKM = "LA6684-0"  # Fragments per kilobase per million (FPKM)
    
    # Systems
    LOINC_SYSTEM = "http://loinc.org"
    HGNC_SYSTEM = "http://www.genenames.org"
    UCUM_SYSTEM = "http://unitsofmeasure.org"
    
    def __init__(self):
        """Initialize FHIR genomics observation creator."""
        logger.info("Initialized FHIR Genomics Observation creator")
    
    def create_gene_expression_observation(
        self,
        gene_id: str,
        gene_symbol: str,
        expression_value: float,
        expression_unit: str = "TPM",
        patient_id: Optional[str] = None,
        specimen_id: Optional[str] = None,
        effective_datetime: Optional[str] = None,
        status: str = "final"
    ) -> Dict[str, Any]:
        """
        Create a FHIR Observation resource for gene expression.
        
        Args:
            gene_id: Gene identifier (e.g., ENSG00000141510)
            gene_symbol: HGNC gene symbol (e.g., TP53)
            expression_value: Expression value (TPM, FPKM, etc.)
            expression_unit: Unit of measurement (TPM, FPKM, counts)
            patient_id: Patient reference ID
            specimen_id: Specimen reference ID
            effective_datetime: Observation date/time (ISO format)
            status: Observation status (final, preliminary, etc.)
            
        Returns:
            FHIR Observation resource as dictionary
        """
        logger.info(f"Creating FHIR observation for gene: {gene_symbol}")
        
        if effective_datetime is None:
            effective_datetime = datetime.now().isoformat()
        
        # Map units to LOINC/UCUM codes
        unit_mapping = {
            "TPM": ("TPM", self.UCUM_SYSTEM, "TPM"),
            "FPKM": ("FPKM", self.UCUM_SYSTEM, "FPKM"),
            "counts": ("counts", self.UCUM_SYSTEM, "{counts}"),
            "log2": ("log2", self.UCUM_SYSTEM, "log2")
        }
        
        unit_display, unit_system, unit_code = unit_mapping.get(
            expression_unit,
            (expression_unit, self.UCUM_SYSTEM, expression_unit)
        )
        
        observation = {
            "resourceType": "Observation",
            "id": f"gene-expression-{gene_id}-{int(datetime.now().timestamp())}",
            "status": status,
            "category": [
                {
                    "coding": [
                        {
                            "system": "http://terminology.hl7.org/CodeSystem/observation-category",
                            "code": "laboratory",
                            "display": "Laboratory"
                        }
                    ]
                }
            ],
            "code": {
                "coding": [
                    {
                        "system": self.LOINC_SYSTEM,
                        "code": self.LOINC_GENE_EXPRESSION,
                        "display": "Gene expression"
                    }
                ],
                "text": f"{gene_symbol} gene expression"
            },
            "effectiveDateTime": effective_datetime,
            "valueQuantity": {
                "value": expression_value,
                "unit": unit_display,
                "system": unit_system,
                "code": unit_code
            },
            "component": [
                {
                    "code": {
                        "coding": [
                            {
                                "system": self.LOINC_SYSTEM,
                                "code": self.LOINC_GENE_ID,
                                "display": "Gene studied [ID]"
                            }
                        ]
                    },
                    "valueCodeableConcept": {
                        "coding": [
                            {
                                "system": self.HGNC_SYSTEM,
                                "code": gene_symbol,
                                "display": gene_symbol
                            }
                        ],
                        "text": gene_id
                    }
                }
            ]
        }
        
        # Add subject (patient) reference if provided
        if patient_id:
            observation["subject"] = {
                "reference": f"Patient/{patient_id}",
                "display": f"Patient {patient_id}"
            }
        
        # Add specimen reference if provided
        if specimen_id:
            observation["specimen"] = {
                "reference": f"Specimen/{specimen_id}",
                "display": f"Specimen {specimen_id}"
            }
        
        return observation
    
    def create_differential_expression_observation(
        self,
        gene_id: str,
        gene_symbol: str,
        log2_fold_change: float,
        p_value: float,
        adjusted_p_value: float,
        base_mean: float,
        patient_id: Optional[str] = None,
        effective_datetime: Optional[str] = None
    ) -> Dict[str, Any]:
        """
        Create FHIR Observation for differential expression analysis results.
        
        Args:
            gene_id: Gene identifier
            gene_symbol: HGNC gene symbol
            log2_fold_change: Log2 fold change value
            p_value: Statistical p-value
            adjusted_p_value: FDR-adjusted p-value
            base_mean: Base mean expression
            patient_id: Patient reference ID
            effective_datetime: Observation date/time
            
        Returns:
            FHIR Observation resource as dictionary
        """
        logger.info(f"Creating differential expression observation for: {gene_symbol}")
        
        if effective_datetime is None:
            effective_datetime = datetime.now().isoformat()
        
        observation = {
            "resourceType": "Observation",
            "id": f"de-gene-{gene_id}-{int(datetime.now().timestamp())}",
            "status": "final",
            "category": [
                {
                    "coding": [
                        {
                            "system": "http://terminology.hl7.org/CodeSystem/observation-category",
                            "code": "laboratory",
                            "display": "Laboratory"
                        }
                    ]
                }
            ],
            "code": {
                "coding": [
                    {
                        "system": self.LOINC_SYSTEM,
                        "code": "81265-1",
                        "display": "Differential gene expression"
                    }
                ],
                "text": f"{gene_symbol} differential expression"
            },
            "effectiveDateTime": effective_datetime,
            "component": [
                {
                    "code": {
                        "coding": [
                            {
                                "system": self.LOINC_SYSTEM,
                                "code": self.LOINC_GENE_ID,
                                "display": "Gene studied [ID]"
                            }
                        ]
                    },
                    "valueCodeableConcept": {
                        "coding": [
                            {
                                "system": self.HGNC_SYSTEM,
                                "code": gene_symbol,
                                "display": gene_symbol
                            }
                        ],
                        "text": gene_id
                    }
                },
                {
                    "code": {
                        "text": "Log2 Fold Change"
                    },
                    "valueQuantity": {
                        "value": log2_fold_change,
                        "unit": "log2",
                        "system": self.UCUM_SYSTEM,
                        "code": "log2"
                    }
                },
                {
                    "code": {
                        "text": "P-value"
                    },
                    "valueQuantity": {
                        "value": p_value,
                        "system": self.UCUM_SYSTEM,
                        "code": "1"
                    }
                },
                {
                    "code": {
                        "text": "Adjusted P-value"
                    },
                    "valueQuantity": {
                        "value": adjusted_p_value,
                        "system": self.UCUM_SYSTEM,
                        "code": "1"
                    }
                },
                {
                    "code": {
                        "text": "Base Mean Expression"
                    },
                    "valueQuantity": {
                        "value": base_mean,
                        "unit": "counts",
                        "system": self.UCUM_SYSTEM,
                        "code": "{counts}"
                    }
                }
            ]
        }
        
        if patient_id:
            observation["subject"] = {
                "reference": f"Patient/{patient_id}",
                "display": f"Patient {patient_id}"
            }
        
        return observation
    
    def serialize_to_json(self, observation: Dict[str, Any], pretty: bool = True) -> str:
        """
        Serialize FHIR observation to JSON string.
        
        Args:
            observation: FHIR observation dictionary
            pretty: Whether to pretty-print JSON
            
        Returns:
            JSON string
        """
        if pretty:
            return json.dumps(observation, indent=2)
        return json.dumps(observation)
    
    def validate_observation(self, observation: Dict[str, Any]) -> bool:
        """
        Basic validation of FHIR observation structure.
        
        Args:
            observation: FHIR observation dictionary
            
        Returns:
            True if valid, False otherwise
        """
        required_fields = ["resourceType", "status", "code"]
        
        for field in required_fields:
            if field not in observation:
                logger.error(f"Missing required field: {field}")
                return False
        
        if observation["resourceType"] != "Observation":
            logger.error("Invalid resourceType")
            return False
        
        logger.info("Observation validation passed")
        return True


# Example usage
if __name__ == "__main__":
    # Create FHIR genomics observation creator
    fhir_creator = FHIRGenomicsObservation()
    
    # Example 1: Gene expression observation
    expression_obs = fhir_creator.create_gene_expression_observation(
        gene_id="ENSG00000141510",
        gene_symbol="TP53",
        expression_value=1234.56,
        expression_unit="TPM",
        patient_id="patient-123",
        specimen_id="specimen-456"
    )
    
    print("=== Gene Expression Observation ===")
    print(fhir_creator.serialize_to_json(expression_obs))
    print()
    
    # Example 2: Differential expression observation
    de_obs = fhir_creator.create_differential_expression_observation(
        gene_id="ENSG00000141510",
        gene_symbol="TP53",
        log2_fold_change=2.5,
        p_value=0.001,
        adjusted_p_value=0.01,
        base_mean=5678.9,
        patient_id="patient-123"
    )
    
    print("=== Differential Expression Observation ===")
    print(fhir_creator.serialize_to_json(de_obs))
    print()
    
    # Validate observations
    print("=== Validation Results ===")
    print(f"Expression observation valid: {fhir_creator.validate_observation(expression_obs)}")
    print(f"DE observation valid: {fhir_creator.validate_observation(de_obs)}")
