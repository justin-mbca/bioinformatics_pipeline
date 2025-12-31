# Healthcare Data Standards Integration

## Overview

This directory contains production-ready implementations for integrating bioinformatics data with standard healthcare data formats used in clinical and research settings.

## Supported Standards

### 1. FHIR R4 (Fast Healthcare Interoperability Resources)

**File:** `fhir_genomics_integration.py`

FHIR is the modern standard for healthcare data exchange, developed by HL7. Our implementation creates FHIR R4 Observation resources for genomic data.

**Features:**
- Gene expression observations with LOINC codes
- Differential expression results
- Patient and specimen references
- HGNC gene identifiers
- Support for TPM, FPKM, and count units
- JSON serialization
- Basic validation

**Use Cases:**
- Integration with EHR systems
- Research data repositories (dbGaP, European Genome-phenome Archive)
- Clinical genomics workflows
- Precision medicine applications

**Example:**
```python
from fhir_genomics_integration import FHIRGenomicsObservation

fhir = FHIRGenomicsObservation()
obs = fhir.create_gene_expression_observation(
    gene_id="ENSG00000141510",
    gene_symbol="TP53",
    expression_value=1234.56,
    expression_unit="TPM",
    patient_id="patient-123"
)
print(fhir.serialize_to_json(obs))
```

### 2. HL7 v2.x Message Parsing

**File:** `hl7_v2_parser.py`

HL7 v2 is the legacy standard still widely used in healthcare systems for laboratory results and patient data.

**Features:**
- Parse ORU^R01 (Lab Results) messages
- Extract patient demographics (PID segment)
- Parse observations (OBX segments)
- Map to RNA-Seq metadata
- Convert to structured JSON
- Error handling for malformed messages

**Use Cases:**
- Laboratory information system (LIS) integration
- Clinical data extraction
- Legacy system integration
- Real-time data feeds

**Example:**
```python
from hl7_v2_parser import HL7v2Parser

parser = HL7v2Parser()
message = parser.parse_message(hl7_text)
metadata = parser.map_to_rnaseq_metadata(message)
```

### 3. OMOP Common Data Model (CDM)

**File:** `omop_cdm_mapping.py`

OMOP CDM is the standard for observational health research, enabling large-scale cohort studies across institutions.

**Features:**
- Map to PERSON table (demographics)
- Map to SPECIMEN table (tissue samples)
- Map to MEASUREMENT table (gene expression)
- Map to OBSERVATION table (clinical metadata)
- Generate SQL INSERT statements
- Concept ID mappings

**Use Cases:**
- Multi-site research studies
- Cohort building
- Comparative effectiveness research
- Real-world evidence generation

**Example:**
```python
from omop_cdm_mapping import OMOPCDMMapper

mapper = OMOPCDMMapper()
person = mapper.map_patient_to_person(
    patient_id="PAT123",
    year_of_birth=1980,
    gender="M"
)
sql = mapper.generate_sql_insert("PERSON", person)
```

## Integration Architecture

```
┌─────────────────┐
│  EHR System     │
│  (HL7 v2/FHIR)  │
└────────┬────────┘
         │
         ▼
┌─────────────────┐      ┌──────────────────┐
│  Data Parser    │─────▶│  RNA-Seq Pipeline│
│  (HL7/FHIR)     │      │  (AWS Glue/EMR)  │
└─────────────────┘      └────────┬─────────┘
                                  │
                                  ▼
                         ┌──────────────────┐
                         │  Results Export  │
                         │  (FHIR/OMOP)     │
                         └────────┬─────────┘
                                  │
                                  ▼
                         ┌──────────────────┐
                         │  Research DB     │
                         │  (OMOP CDM)      │
                         └──────────────────┘
```

## Concept ID Mappings

### OMOP Standard Concepts

| Concept | Concept ID | Description |
|---------|-----------|-------------|
| Male | 8507 | Gender: Male |
| Female | 8532 | Gender: Female |
| Blood specimen | 32856 | Specimen type |
| Tissue specimen | 32849 | Specimen type |
| Lab test | 44818702 | Measurement type |

### Custom Genomics Concepts

| Concept | Concept ID | Description |
|---------|-----------|-------------|
| Gene expression TPM | 2000000001 | TPM measurement |
| Gene expression FPKM | 2000000002 | FPKM measurement |
| Gene expression count | 2000000003 | Raw count |
| Log2 fold change | 2000000004 | Differential expression |

**Note:** Custom concept IDs should be registered in your OMOP vocabulary.

## LOINC Codes for Genomics

| LOINC Code | Description |
|-----------|-------------|
| 81247-9 | Gene expression |
| 48018-6 | Gene studied [ID] |
| 81265-1 | Differential gene expression |
| LA6683-2 | Transcripts per million (TPM) |
| LA6684-0 | Fragments per kilobase per million (FPKM) |

## Validation and Testing

Run validation scripts:

```bash
# Test FHIR integration
python fhir_genomics_integration.py

# Test HL7 parsing
python hl7_v2_parser.py

# Test OMOP mapping
python omop_cdm_mapping.py
```

## Compliance Considerations

### HIPAA Compliance

- **De-identification:** Remove or encrypt PHI/PII before processing
- **Access Controls:** Implement role-based access
- **Audit Logging:** Track all data access
- **Encryption:** Encrypt data at rest and in transit
- **Business Associate Agreements:** Ensure vendor compliance

### FDA 21 CFR Part 11

- **Electronic Signatures:** Implement for data approval
- **Audit Trails:** Immutable records of data changes
- **System Validation:** Document software validation
- **Data Integrity:** Checksums and version control

### GDPR Compliance

- **Right to Access:** Provide data exports
- **Right to Erasure:** Implement data deletion
- **Data Portability:** Support standard formats (FHIR)
- **Privacy by Design:** Minimize data collection

## Best Practices

1. **Data Quality:**
   - Validate all incoming data
   - Handle missing values appropriately
   - Document data transformations

2. **Interoperability:**
   - Use standard vocabularies (LOINC, SNOMED CT, HGNC)
   - Map to standard concept IDs
   - Version your data schemas

3. **Security:**
   - Never log PHI/PII
   - Encrypt sensitive data
   - Implement access controls
   - Regular security audits

4. **Documentation:**
   - Document all mappings
   - Maintain data dictionaries
   - Version control configurations

## Common Integration Patterns

### Pattern 1: EHR to RNA-Seq Pipeline

```python
# 1. Receive HL7 message from EHR
hl7_message = receive_hl7_message()

# 2. Parse and extract metadata
parser = HL7v2Parser()
parsed = parser.parse_message(hl7_message)
metadata = parser.map_to_rnaseq_metadata(parsed)

# 3. Trigger RNA-Seq pipeline with metadata
trigger_glue_job(metadata)
```

### Pattern 2: Results to Research Database

```python
# 1. Get DE results
de_results = read_de_results_from_s3()

# 2. Map to OMOP CDM
mapper = OMOPCDMMapper()
for gene in de_results:
    measurement = mapper.map_gene_expression_to_measurement(
        person_id=patient_omop_id,
        gene_symbol=gene['symbol'],
        expression_value=gene['tpm'],
        measurement_date=date.today()
    )
    sql = mapper.generate_sql_insert("MEASUREMENT", measurement)
    execute_sql(sql)
```

### Pattern 3: FHIR API Export

```python
# 1. Get significant genes
sig_genes = get_significant_genes()

# 2. Create FHIR observations
fhir = FHIRGenomicsObservation()
observations = []
for gene in sig_genes:
    obs = fhir.create_differential_expression_observation(
        gene_id=gene['gene_id'],
        gene_symbol=gene['symbol'],
        log2_fold_change=gene['log2fc'],
        p_value=gene['pvalue'],
        adjusted_p_value=gene['padj'],
        base_mean=gene['baseMean']
    )
    observations.append(obs)

# 3. POST to FHIR server
for obs in observations:
    post_to_fhir_server(obs)
```

## Troubleshooting

### Common Issues

**HL7 Parsing Errors:**
- Check field separator (usually |)
- Verify segment order (MSH must be first)
- Handle escape characters properly

**FHIR Validation Errors:**
- Ensure required fields are present
- Validate against FHIR schema
- Check concept code systems

**OMOP Concept ID Issues:**
- Verify concept IDs exist in vocabulary
- Use standard concepts when possible
- Document custom concept mappings

## Additional Resources

- [FHIR Genomics Implementation Guide](http://hl7.org/fhir/uv/genomics-reporting/)
- [HL7 v2 Specification](http://www.hl7.org/implement/standards/product_brief.cfm?product_id=185)
- [OMOP CDM Documentation](https://ohdsi.github.io/CommonDataModel/)
- [LOINC Database](https://loinc.org/)
- [HGNC Database](https://www.genenames.org/)
- [HIPAA Guidelines](https://www.hhs.gov/hipaa/)
- [FDA 21 CFR Part 11](https://www.fda.gov/regulatory-information/search-fda-guidance-documents/part-11-electronic-records-electronic-signatures-scope-and-application)

## Support and Contact

For questions or issues with healthcare data standard integration:
- Open an issue in the repository
- Contact the bioinformatics team
- Refer to the documentation

## License

This code is provided for educational and research purposes. Consult with legal and compliance teams before using in production healthcare settings.
