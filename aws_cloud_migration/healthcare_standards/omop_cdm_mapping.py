"""
OMOP Common Data Model (CDM) Mapping

This module maps bioinformatics data (gene expression, patient demographics,
specimen information) to the OMOP CDM v5.4 standard tables for research cohort analysis.
"""

import logging
from typing import Dict, List, Optional, Any
from datetime import datetime, date
from dataclasses import dataclass

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@dataclass
class OMOPPerson:
    """OMOP PERSON table record."""
    person_id: int
    gender_concept_id: int
    year_of_birth: int
    month_of_birth: Optional[int] = None
    day_of_birth: Optional[int] = None
    race_concept_id: int = 0
    ethnicity_concept_id: int = 0
    location_id: Optional[int] = None
    provider_id: Optional[int] = None
    care_site_id: Optional[int] = None
    person_source_value: str = ""
    gender_source_value: str = ""


@dataclass
class OMOPSpecimen:
    """OMOP SPECIMEN table record."""
    specimen_id: int
    person_id: int
    specimen_concept_id: int
    specimen_type_concept_id: int
    specimen_date: date
    specimen_datetime: Optional[datetime] = None
    quantity: Optional[float] = None
    unit_concept_id: Optional[int] = None
    anatomic_site_concept_id: Optional[int] = None
    disease_status_concept_id: Optional[int] = None
    specimen_source_id: str = ""
    specimen_source_value: str = ""
    unit_source_value: Optional[str] = None
    anatomic_site_source_value: Optional[str] = None
    disease_status_source_value: Optional[str] = None


@dataclass
class OMOPMeasurement:
    """OMOP MEASUREMENT table record."""
    measurement_id: int
    person_id: int
    measurement_concept_id: int
    measurement_date: date
    measurement_datetime: Optional[datetime] = None
    measurement_time: Optional[str] = None
    measurement_type_concept_id: int = 0
    operator_concept_id: Optional[int] = None
    value_as_number: Optional[float] = None
    value_as_concept_id: Optional[int] = None
    unit_concept_id: Optional[int] = None
    range_low: Optional[float] = None
    range_high: Optional[float] = None
    provider_id: Optional[int] = None
    visit_occurrence_id: Optional[int] = None
    visit_detail_id: Optional[int] = None
    measurement_source_value: str = ""
    measurement_source_concept_id: Optional[int] = None
    unit_source_value: Optional[str] = None
    value_source_value: Optional[str] = None


@dataclass
class OMOPObservation:
    """OMOP OBSERVATION table record."""
    observation_id: int
    person_id: int
    observation_concept_id: int
    observation_date: date
    observation_datetime: Optional[datetime] = None
    observation_type_concept_id: int = 0
    value_as_number: Optional[float] = None
    value_as_string: Optional[str] = None
    value_as_concept_id: Optional[int] = None
    qualifier_concept_id: Optional[int] = None
    unit_concept_id: Optional[int] = None
    provider_id: Optional[int] = None
    visit_occurrence_id: Optional[int] = None
    visit_detail_id: Optional[int] = None
    observation_source_value: str = ""
    observation_source_concept_id: Optional[int] = None
    unit_source_value: Optional[str] = None
    qualifier_source_value: Optional[str] = None


class OMOPCDMMapper:
    """
    Maps bioinformatics data to OMOP CDM tables.
    
    Key OMOP concept IDs for genomics:
    - Gender: Male=8507, Female=8532, Unknown=0
    - Specimen Type: Blood=32856, Tissue=32849
    - Measurement Type: Lab Test=44818702
    - Unit: TPM (custom), FPKM (custom), counts (custom)
    """
    
    # Standard OMOP Concept IDs
    CONCEPT_GENDER_MALE = 8507
    CONCEPT_GENDER_FEMALE = 8532
    CONCEPT_GENDER_UNKNOWN = 0
    
    CONCEPT_SPECIMEN_BLOOD = 32856
    CONCEPT_SPECIMEN_TISSUE = 32849
    CONCEPT_SPECIMEN_TYPE_BIOBANK = 32850
    
    CONCEPT_MEASUREMENT_TYPE_LAB = 44818702
    CONCEPT_OBSERVATION_TYPE_EHR = 38000280
    
    # Custom concept IDs for genomics (would be in custom vocabulary)
    CONCEPT_GENE_EXPRESSION_TPM = 2000000001
    CONCEPT_GENE_EXPRESSION_FPKM = 2000000002
    CONCEPT_GENE_EXPRESSION_COUNT = 2000000003
    CONCEPT_LOG2_FOLD_CHANGE = 2000000004
    
    def __init__(self):
        """Initialize OMOP CDM mapper."""
        logger.info("Initialized OMOP CDM mapper")
        self.person_counter = 1
        self.specimen_counter = 1
        self.measurement_counter = 1
        self.observation_counter = 1
    
    def map_patient_to_person(
        self,
        patient_id: str,
        year_of_birth: int,
        gender: str,
        month_of_birth: Optional[int] = None,
        day_of_birth: Optional[int] = None
    ) -> OMOPPerson:
        """
        Map patient demographics to OMOP PERSON table.
        
        Args:
            patient_id: Source patient identifier
            year_of_birth: Year of birth
            gender: Gender (M/F/U)
            month_of_birth: Month of birth (1-12)
            day_of_birth: Day of birth (1-31)
            
        Returns:
            OMOPPerson record
        """
        logger.info(f"Mapping patient {patient_id} to OMOP PERSON")
        
        # Map gender to concept ID
        gender_concept_map = {
            'M': self.CONCEPT_GENDER_MALE,
            'F': self.CONCEPT_GENDER_FEMALE,
            'U': self.CONCEPT_GENDER_UNKNOWN
        }
        gender_concept_id = gender_concept_map.get(gender.upper(), self.CONCEPT_GENDER_UNKNOWN)
        
        person = OMOPPerson(
            person_id=self.person_counter,
            gender_concept_id=gender_concept_id,
            year_of_birth=year_of_birth,
            month_of_birth=month_of_birth,
            day_of_birth=day_of_birth,
            person_source_value=patient_id,
            gender_source_value=gender
        )
        
        self.person_counter += 1
        return person
    
    def map_tissue_sample_to_specimen(
        self,
        specimen_id: str,
        person_id: int,
        collection_date: date,
        tissue_type: str = "tissue",
        anatomic_site: Optional[str] = None
    ) -> OMOPSpecimen:
        """
        Map tissue sample to OMOP SPECIMEN table.
        
        Args:
            specimen_id: Source specimen identifier
            person_id: OMOP person_id
            collection_date: Date of collection
            tissue_type: Type of specimen (blood, tissue, etc.)
            anatomic_site: Anatomical site description
            
        Returns:
            OMOPSpecimen record
        """
        logger.info(f"Mapping specimen {specimen_id} to OMOP SPECIMEN")
        
        # Map specimen type
        specimen_concept_map = {
            'blood': self.CONCEPT_SPECIMEN_BLOOD,
            'tissue': self.CONCEPT_SPECIMEN_TISSUE
        }
        specimen_concept_id = specimen_concept_map.get(tissue_type.lower(), self.CONCEPT_SPECIMEN_TISSUE)
        
        specimen = OMOPSpecimen(
            specimen_id=self.specimen_counter,
            person_id=person_id,
            specimen_concept_id=specimen_concept_id,
            specimen_type_concept_id=self.CONCEPT_SPECIMEN_TYPE_BIOBANK,
            specimen_date=collection_date,
            specimen_datetime=datetime.combine(collection_date, datetime.min.time()),
            specimen_source_id=specimen_id,
            specimen_source_value=tissue_type,
            anatomic_site_source_value=anatomic_site
        )
        
        self.specimen_counter += 1
        return specimen
    
    def map_gene_expression_to_measurement(
        self,
        person_id: int,
        gene_symbol: str,
        expression_value: float,
        measurement_date: date,
        expression_unit: str = "TPM",
        specimen_id: Optional[int] = None
    ) -> OMOPMeasurement:
        """
        Map gene expression value to OMOP MEASUREMENT table.
        
        Args:
            person_id: OMOP person_id
            gene_symbol: Gene symbol (e.g., TP53)
            expression_value: Expression value
            measurement_date: Date of measurement
            expression_unit: Unit (TPM, FPKM, counts)
            specimen_id: Optional specimen_id
            
        Returns:
            OMOPMeasurement record
        """
        logger.info(f"Mapping gene expression {gene_symbol} to OMOP MEASUREMENT")
        
        # Map unit to concept ID
        unit_concept_map = {
            'TPM': self.CONCEPT_GENE_EXPRESSION_TPM,
            'FPKM': self.CONCEPT_GENE_EXPRESSION_FPKM,
            'counts': self.CONCEPT_GENE_EXPRESSION_COUNT
        }
        unit_concept_id = unit_concept_map.get(expression_unit, self.CONCEPT_GENE_EXPRESSION_TPM)
        
        measurement = OMOPMeasurement(
            measurement_id=self.measurement_counter,
            person_id=person_id,
            measurement_concept_id=unit_concept_id,
            measurement_date=measurement_date,
            measurement_datetime=datetime.combine(measurement_date, datetime.min.time()),
            measurement_type_concept_id=self.CONCEPT_MEASUREMENT_TYPE_LAB,
            value_as_number=expression_value,
            unit_concept_id=unit_concept_id,
            measurement_source_value=f"{gene_symbol}_expression",
            unit_source_value=expression_unit
        )
        
        self.measurement_counter += 1
        return measurement
    
    def map_differential_expression_to_observation(
        self,
        person_id: int,
        gene_symbol: str,
        log2_fold_change: float,
        p_value: float,
        adjusted_p_value: float,
        observation_date: date
    ) -> OMOPObservation:
        """
        Map differential expression result to OMOP OBSERVATION table.
        
        Args:
            person_id: OMOP person_id
            gene_symbol: Gene symbol
            log2_fold_change: Log2 fold change value
            p_value: Statistical p-value
            adjusted_p_value: FDR-adjusted p-value
            observation_date: Date of observation
            
        Returns:
            OMOPObservation record
        """
        logger.info(f"Mapping differential expression {gene_symbol} to OMOP OBSERVATION")
        
        # Create observation with log2FC as value
        observation = OMOPObservation(
            observation_id=self.observation_counter,
            person_id=person_id,
            observation_concept_id=self.CONCEPT_LOG2_FOLD_CHANGE,
            observation_date=observation_date,
            observation_datetime=datetime.combine(observation_date, datetime.min.time()),
            observation_type_concept_id=self.CONCEPT_OBSERVATION_TYPE_EHR,
            value_as_number=log2_fold_change,
            value_as_string=f"gene={gene_symbol};pval={p_value};padj={adjusted_p_value}",
            observation_source_value=f"{gene_symbol}_log2fc"
        )
        
        self.observation_counter += 1
        return observation
    
    def generate_sql_insert(self, table_name: str, record: Any) -> str:
        """
        Generate SQL INSERT statement for OMOP CDM table.
        
        Args:
            table_name: OMOP table name
            record: Record dataclass instance
            
        Returns:
            SQL INSERT statement
        """
        # Convert dataclass to dict
        if hasattr(record, '__dataclass_fields__'):
            data = {k: v for k, v in record.__dict__.items() if v is not None}
        else:
            data = record
        
        columns = ', '.join(data.keys())
        
        # Format values
        values = []
        for v in data.values():
            if isinstance(v, str):
                values.append(f"'{v}'")
            elif isinstance(v, (date, datetime)):
                values.append(f"'{v.isoformat()}'")
            elif v is None:
                values.append('NULL')
            else:
                values.append(str(v))
        
        values_str = ', '.join(values)
        
        sql = f"INSERT INTO {table_name} ({columns}) VALUES ({values_str});"
        return sql


# Example usage
if __name__ == "__main__":
    mapper = OMOPCDMMapper()
    
    # Example 1: Map patient demographics
    print("=== OMOP PERSON Mapping ===")
    person = mapper.map_patient_to_person(
        patient_id="PAT123456",
        year_of_birth=1980,
        gender="M",
        month_of_birth=1,
        day_of_birth=15
    )
    print(mapper.generate_sql_insert("PERSON", person))
    print()
    
    # Example 2: Map specimen
    print("=== OMOP SPECIMEN Mapping ===")
    specimen = mapper.map_tissue_sample_to_specimen(
        specimen_id="SPEC789",
        person_id=person.person_id,
        collection_date=date(2024, 1, 1),
        tissue_type="tissue",
        anatomic_site="lung"
    )
    print(mapper.generate_sql_insert("SPECIMEN", specimen))
    print()
    
    # Example 3: Map gene expression
    print("=== OMOP MEASUREMENT Mapping ===")
    measurement = mapper.map_gene_expression_to_measurement(
        person_id=person.person_id,
        gene_symbol="TP53",
        expression_value=1234.56,
        measurement_date=date(2024, 1, 1),
        expression_unit="TPM"
    )
    print(mapper.generate_sql_insert("MEASUREMENT", measurement))
    print()
    
    # Example 4: Map differential expression
    print("=== OMOP OBSERVATION Mapping ===")
    observation = mapper.map_differential_expression_to_observation(
        person_id=person.person_id,
        gene_symbol="TP53",
        log2_fold_change=2.5,
        p_value=0.001,
        adjusted_p_value=0.01,
        observation_date=date(2024, 1, 1)
    )
    print(mapper.generate_sql_insert("OBSERVATION", observation))
    print()
    
    # Example 5: Document concept mappings
    print("=== Concept ID Mappings ===")
    print(f"Gender Male: {mapper.CONCEPT_GENDER_MALE}")
    print(f"Gender Female: {mapper.CONCEPT_GENDER_FEMALE}")
    print(f"Specimen Tissue: {mapper.CONCEPT_SPECIMEN_TISSUE}")
    print(f"Gene Expression TPM: {mapper.CONCEPT_GENE_EXPRESSION_TPM}")
    print(f"Log2 Fold Change: {mapper.CONCEPT_LOG2_FOLD_CHANGE}")
