# AWS S3 Data Lake for Bioinformatics Pipeline

## Overview

This directory contains configurations and scripts for managing the S3 Data Lake that serves as the central storage layer for the bioinformatics pipeline.

## S3 Bucket Structure

```
s3://bioinformatics-pipeline/
├── raw/
│   ├── counts/                    # Raw RNA-Seq count matrices (CSV/TSV)
│   │   ├── study01/
│   │   ├── study02/
│   │   └── archived/
│   ├── metadata/                  # Sample metadata and clinical data
│   │   ├── sample_info.csv
│   │   ├── patient_data.csv
│   │   └── experimental_design.json
│   └── fastq/                     # Raw sequencing reads (optional)
│       ├── paired_end/
│       └── single_end/
├── processed/
│   ├── normalized/                # Normalized count data (Parquet)
│   │   ├── date=2024-01-01/
│   │   └── date=2024-01-02/
│   ├── filtered/                  # QC-filtered data
│   └── transformed/               # Additional transformations
├── results/
│   ├── de_genes/                  # Differential expression results
│   │   ├── treatment_vs_control/
│   │   └── timepoint_comparisons/
│   ├── enrichment/                # Pathway enrichment analysis
│   ├── clustering/                # Clustering results
│   └── ml_models/                 # Machine learning models
├── logs/
│   ├── glue/                      # AWS Glue job logs
│   ├── emr/                       # EMR cluster logs
│   ├── lambda/                    # Lambda function logs
│   └── ecs/                       # ECS task logs
├── catalog/                       # Glue Data Catalog metadata
├── glue-scripts/                  # Glue ETL scripts
├── glue-temp/                     # Temporary files for Glue jobs
├── bootstrap/                     # EMR bootstrap scripts
└── emr-scripts/                   # EMR Spark scripts
```

## Data Organization Best Practices

### 1. Naming Conventions

**File Naming:**
- Use lowercase with underscores
- Include date stamps: `YYYYMMDD` or `YYYY-MM-DD`
- Include version if applicable: `v1`, `v2`
- Be descriptive: `study01_rnaseq_counts_20240101.csv`

**Directory Naming:**
- Use descriptive, hierarchical names
- Separate by project, experiment, or date
- Use hyphens or underscores, not spaces

**Examples:**
```
raw/counts/lung_cancer_study/cohort_a_counts_20240101.csv
processed/normalized/date=2024-01-01/lung_cancer_normalized.parquet
results/de_genes/treated_vs_control_20240101/significant_genes.parquet
```

### 2. Partitioning Strategy

**Time-based Partitioning:**
```
processed/normalized/
├── date=2024-01-01/
├── date=2024-01-02/
└── date=2024-01-03/
```

**Study-based Partitioning:**
```
raw/counts/
├── study=lung_cancer/
├── study=breast_cancer/
└── study=melanoma/
```

**Multi-level Partitioning:**
```
results/de_genes/
├── study=lung_cancer/
│   └── comparison=treated_vs_control/
│       └── date=2024-01-01/
└── study=breast_cancer/
    └── comparison=mutant_vs_wildtype/
        └── date=2024-01-01/
```

### 3. File Formats

**Raw Data:**
- CSV/TSV: For count matrices and metadata
- FASTQ: For raw sequencing reads (gzipped)

**Processed Data:**
- Parquet: Primary format for processed data
  - Columnar format for efficient querying
  - Compression support (Snappy, Gzip)
  - Schema evolution support

**Results:**
- Parquet: For analysis results
- JSON: For metadata and configurations
- PDF/PNG: For plots and reports

### 4. Data Lifecycle

**Hot Data (Standard Storage):**
- Recently uploaded raw data (< 30 days)
- Active analysis results
- Frequently accessed datasets

**Warm Data (Standard-IA):**
- Older raw data (30-90 days)
- Completed analysis results
- Reference datasets

**Cold Data (Glacier):**
- Archived raw data (> 90 days)
- Historical analysis results
- Long-term backup data

**Logs (Automatic Deletion):**
- Application logs (90 days retention)
- Temporary files (7 days retention)

## Data Tagging Strategy

Apply consistent tags to S3 objects for cost tracking and data management:

```python
tags = {
    'Project': 'bioinformatics-pipeline',
    'DataType': 'rnaseq',           # rnaseq, scrnaseq, wgs, etc.
    'Classification': 'sensitive',   # public, internal, sensitive, confidential
    'Owner': 'data-science-team',
    'Environment': 'production',     # dev, staging, production
    'CostCenter': 'research-dev',
    'RetentionPeriod': '7years',     # Compliance requirement
    'Study': 'lung-cancer-2024'
}
```

## Access Patterns

### Common Queries

**Athena Query - Recent Results:**
```sql
SELECT gene_name, log2FoldChange, padj
FROM rnaseq_db.de_results
WHERE date >= '2024-01-01'
  AND padj < 0.05
ORDER BY ABS(log2FoldChange) DESC
LIMIT 100;
```

**Spark Query - Aggregate Analysis:**
```python
df = spark.read.parquet("s3://bioinformatics-pipeline/processed/normalized/")
df.filter(df.baseMean > 10).groupBy("chromosome").count().show()
```

### Optimization Tips

1. **Partition Pruning:** Query only necessary partitions
2. **Column Selection:** Select only required columns
3. **Predicate Pushdown:** Apply filters early
4. **Compression:** Use Snappy for balanced performance
5. **File Size:** Target 128-512 MB per file

## Security Considerations

### 1. Encryption

**At Rest:**
- S3 bucket encryption: AES-256 or AWS KMS
- Enable default encryption for bucket

**In Transit:**
- TLS/SSL for all S3 API calls
- HTTPS endpoints only

### 2. Access Control

**Bucket Policy:**
```json
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Sid": "DenyInsecureTransport",
      "Effect": "Deny",
      "Principal": "*",
      "Action": "s3:*",
      "Resource": [
        "arn:aws:s3:::bioinformatics-pipeline",
        "arn:aws:s3:::bioinformatics-pipeline/*"
      ],
      "Condition": {
        "Bool": {
          "aws:SecureTransport": "false"
        }
      }
    }
  ]
}
```

**IAM Policies:**
- Least privilege access
- Separate read and write permissions
- Use IAM roles for services

### 3. Data Classification

- **Public:** No restrictions (rare for genomics data)
- **Internal:** Within organization only
- **Sensitive:** Requires approval, audit logging
- **Confidential:** Encrypted, restricted access, PHI/PII

### 4. Compliance

**HIPAA:**
- Enable S3 access logging
- Use KMS encryption with audit trails
- Implement data retention policies
- Restrict public access

**GDPR:**
- Data subject access rights
- Right to erasure (delete data)
- Data minimization
- Privacy by design

## Monitoring and Alerts

### CloudWatch Metrics

Monitor key S3 metrics:
- Bucket size
- Number of objects
- Request count and latency
- 4xx and 5xx errors

### CloudWatch Alarms

Set up alarms for:
- Bucket size exceeds threshold
- High error rates
- Unusual access patterns
- Failed uploads

### S3 Access Logging

Enable logging to track:
- Who accessed what data
- When access occurred
- Source IP addresses
- Request types (GET, PUT, DELETE)

## Cost Optimization

### Storage Costs (us-east-1)

| Storage Class | Price (per GB/month) | Use Case |
|--------------|---------------------|----------|
| Standard | $0.023 | Active data |
| Standard-IA | $0.0125 | Infrequent access |
| Glacier Flexible | $0.004 | Long-term archive |
| Glacier Deep Archive | $0.00099 | Compliance archives |

### Data Transfer Costs

- **In:** Free
- **Out to Internet:** $0.09 per GB (first 10 TB)
- **Out to CloudFront:** Free
- **Out to AWS services (same region):** Free

### Optimization Strategies

1. **Lifecycle Policies:** Automatic tier transitions
2. **Compression:** Reduce storage size by 70-90%
3. **Delete Old Logs:** Remove temporary files
4. **Intelligent-Tiering:** Automatic cost optimization
5. **S3 Analytics:** Monitor access patterns
6. **Request Optimization:** Batch small files

## Usage Examples

See `upload_datasets.py` for programmatic upload examples.

### AWS CLI Upload

```bash
# Upload single file
aws s3 cp local_counts.csv s3://bioinformatics-pipeline/raw/counts/

# Upload directory
aws s3 sync local_data/ s3://bioinformatics-pipeline/raw/counts/ \
    --exclude "*.tmp" \
    --include "*.csv"

# Upload with metadata
aws s3 cp local_counts.csv s3://bioinformatics-pipeline/raw/counts/ \
    --metadata study=lung_cancer,date=2024-01-01 \
    --storage-class STANDARD

# Upload with server-side encryption
aws s3 cp local_counts.csv s3://bioinformatics-pipeline/raw/counts/ \
    --sse AES256
```

### Python Boto3 Upload

```python
import boto3

s3 = boto3.client('s3')

# Upload file
s3.upload_file(
    'local_counts.csv',
    'bioinformatics-pipeline',
    'raw/counts/study01_counts.csv',
    ExtraArgs={
        'Metadata': {
            'study': 'lung_cancer',
            'date': '2024-01-01'
        },
        'ServerSideEncryption': 'AES256',
        'StorageClass': 'STANDARD'
    }
)
```

## Disaster Recovery

### Backup Strategy

1. **Cross-Region Replication:**
   - Enable for critical data
   - Replicate to different AWS region
   - Automatic and continuous

2. **Versioning:**
   - Enable bucket versioning
   - Protect against accidental deletion
   - Recover previous versions

3. **Lifecycle Policies:**
   - Transition old versions to Glacier
   - Expire noncurrent versions after 90 days

### Recovery Procedures

**Accidental Deletion:**
```bash
# List deleted objects
aws s3api list-object-versions \
    --bucket bioinformatics-pipeline \
    --prefix raw/counts/

# Restore specific version
aws s3api copy-object \
    --copy-source bioinformatics-pipeline/raw/counts/file.csv?versionId=<version_id> \
    --bucket bioinformatics-pipeline \
    --key raw/counts/file.csv
```

## Integration Points

- **AWS Glue:** Data catalog and ETL jobs
- **Amazon Athena:** SQL queries on S3 data
- **AWS EMR:** Large-scale Spark processing
- **AWS Lambda:** Event-driven processing
- **Amazon SageMaker:** ML model training
- **Amazon QuickSight:** BI dashboards

## Additional Resources

- [S3 Best Practices](https://docs.aws.amazon.com/AmazonS3/latest/userguide/best-practices.html)
- [S3 Storage Classes](https://aws.amazon.com/s3/storage-classes/)
- [S3 Lifecycle Management](https://docs.aws.amazon.com/AmazonS3/latest/userguide/object-lifecycle-mgmt.html)
- [S3 Security Best Practices](https://docs.aws.amazon.com/AmazonS3/latest/userguide/security-best-practices.html)
