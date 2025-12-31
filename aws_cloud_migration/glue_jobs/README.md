# AWS Glue ETL Jobs

This directory contains AWS Glue ETL jobs for processing RNA-Seq count data in a serverless, scalable manner.

## Overview

AWS Glue is a fully managed ETL (Extract, Transform, Load) service that makes it easy to prepare and load data for analytics. These scripts demonstrate production-ready patterns for genomics data processing.

## Files

### `etl_rnaseq_to_s3.py`
Complete AWS Glue job that:
- Reads RNA-Seq count data from the AWS Glue Data Catalog
- Filters low-quality reads (count > 10 by default)
- Writes processed data to S3 in Parquet format
- Includes comprehensive error handling and logging
- Uses PySpark for distributed processing

**Key Features:**
- Type hints for better code documentation
- Modular function design for testability
- CloudWatch logging integration
- Job bookmarking support for incremental processing

### `glue_job_definition.json`
Production-ready Glue job configuration including:
- Job name and IAM role
- Script location in S3
- Worker type (G.1X) and count (5 workers)
- Timeout settings (48 hours)
- Spark UI and continuous logging enabled
- Job parameters for database, table, and output paths

### `data_catalog_setup.py`
Boto3 script to programmatically create and manage the Glue Data Catalog:
- Creates `rnaseq_db` database
- Defines table schema for raw RNA-Seq counts
- Supports both CSV and Parquet table formats
- Includes methods to list databases and tables

**Schema Definition:**
- gene_id (string): Ensembl gene identifier
- gene_name (string): HGNC gene symbol
- chromosome (string): Chromosome location
- start_position (int): Gene start position
- end_position (int): Gene end position
- strand (string): Strand orientation (+/-)
- gene_length (int): Gene length in bases
- sample_XX (int): Read counts for each sample

## Prerequisites

1. **AWS Account**: Active AWS account with appropriate permissions
2. **IAM Role**: Create an IAM role for Glue with permissions for:
   - S3 read/write access
   - Glue Data Catalog access
   - CloudWatch Logs write access
3. **S3 Bucket**: Create an S3 bucket for data and scripts:
   ```
   s3://bioinformatics-pipeline/
   ├── glue-scripts/
   ├── glue-temp/
   ├── glue-logs/
   ├── raw/counts/
   └── processed/normalized/
   ```

## Setup Instructions

### 1. Create Glue Data Catalog

Run the catalog setup script to create the database and tables:

```bash
python data_catalog_setup.py
```

This will create:
- Database: `rnaseq_db`
- Table: `raw_counts` (CSV format)
- Table: `processed_counts` (Parquet format)

### 2. Upload ETL Script to S3

```bash
aws s3 cp etl_rnaseq_to_s3.py s3://bioinformatics-pipeline/glue-scripts/
```

### 3. Create the Glue Job

Update the placeholders in `glue_job_definition.json`:
- Replace `<ACCOUNT_ID>` with your AWS account ID
- Replace `<BUCKET_NAME>` with your S3 bucket name

Create the job using AWS CLI:

```bash
aws glue create-job --cli-input-json file://glue_job_definition.json
```

Or use the AWS Console:
1. Navigate to AWS Glue → ETL → Jobs
2. Click "Add job"
3. Upload the JSON configuration

### 4. Upload Sample Data

Upload your RNA-Seq count data to the raw location:

```bash
aws s3 cp sample_counts.csv s3://bioinformatics-pipeline/raw/counts/
```

### 5. Run the ETL Job

Via AWS CLI:
```bash
aws glue start-job-run \
    --job-name rnaseq-etl-pipeline \
    --arguments '{"--count_threshold":"15"}'
```

Via AWS Console:
1. Navigate to AWS Glue → ETL → Jobs
2. Select `rnaseq-etl-pipeline`
3. Click "Run job"

## Configuration

### Job Parameters

You can override default parameters when starting the job:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `database_name` | `rnaseq_db` | Glue database name |
| `table_name` | `raw_counts` | Input table name |
| `output_path` | `s3://...` | S3 output location |
| `count_threshold` | `10` | Minimum read count filter |

Example with custom threshold:
```bash
aws glue start-job-run \
    --job-name rnaseq-etl-pipeline \
    --arguments '{
        "--database_name": "rnaseq_db",
        "--table_name": "raw_counts",
        "--output_path": "s3://bioinformatics-pipeline/processed/normalized/",
        "--count_threshold": "20"
    }'
```

## Monitoring

### CloudWatch Logs

View job execution logs in CloudWatch:
1. Navigate to CloudWatch → Log groups
2. Find log group: `/aws-glue/jobs/output`
3. Select the log stream for your job run

### Spark UI

Access the Spark UI for detailed execution metrics:
1. Navigate to AWS Glue → ETL → Jobs
2. Select your job → History tab
3. Click "View Spark UI" for a specific run

### Metrics

Key metrics to monitor:
- **Data Processing Metrics**: Bytes read, bytes written, rows processed
- **Performance Metrics**: Execution time, DPU hours
- **Error Metrics**: Failed tasks, retries

## Troubleshooting

### Job Fails with "Table Not Found"

Ensure the Glue Data Catalog is properly configured:
```bash
aws glue get-table --database-name rnaseq_db --name raw_counts
```

### Out of Memory Errors

Increase worker capacity in job definition:
```json
{
  "WorkerType": "G.2X",
  "NumberOfWorkers": 10
}
```

### S3 Access Denied

Verify IAM role permissions:
```json
{
  "Effect": "Allow",
  "Action": [
    "s3:GetObject",
    "s3:PutObject",
    "s3:ListBucket"
  ],
  "Resource": [
    "arn:aws:s3:::bioinformatics-pipeline/*"
  ]
}
```

## Cost Optimization

- **Job Bookmarks**: Enable to process only new data
- **Worker Type**: Start with G.1X, scale to G.2X if needed
- **Worker Count**: Monitor DPU utilization and adjust
- **Timeout**: Set appropriate timeout to avoid runaway costs
- **Scheduling**: Use triggers for event-driven processing

**Estimated Costs** (us-east-1):
- G.1X worker: $0.44/hour per DPU
- 5 workers = 5 DPU
- Processing 100GB RNA-Seq data: ~2 hours = $4.40

## Security Best Practices

1. **Encryption**: Enable S3 encryption at rest (AES-256 or KMS)
2. **IAM Roles**: Use least-privilege permissions
3. **VPC**: Run Glue jobs in VPC for private data access
4. **Network**: Use VPC endpoints for S3 access
5. **Logging**: Enable CloudTrail for audit logs
6. **Data Classification**: Tag sensitive data appropriately

## Integration with Other Services

### Trigger from S3 Events (via Lambda)
See `../lambda_functions/trigger_pipeline.py` for automated job triggering

### Output to EMR
Processed Parquet files can be consumed by EMR Spark jobs for downstream analysis

### Query with Athena
Query processed data directly using Amazon Athena:
```sql
SELECT gene_name, AVG(sample_01) as avg_count
FROM rnaseq_db.processed_counts
WHERE chromosome = 'chr1'
GROUP BY gene_name
ORDER BY avg_count DESC
LIMIT 10;
```

## Best Practices

1. **Partitioning**: Partition large datasets by date or chromosome
2. **Compression**: Use Snappy compression for Parquet files
3. **Column Pruning**: Select only required columns
4. **Predicate Pushdown**: Filter data early in the pipeline
5. **Testing**: Test with small datasets before full production runs
6. **Documentation**: Document schema changes and processing logic

## Additional Resources

- [AWS Glue Documentation](https://docs.aws.amazon.com/glue/)
- [AWS Glue Best Practices](https://docs.aws.amazon.com/glue/latest/dg/best-practices.html)
- [PySpark API Reference](https://spark.apache.org/docs/latest/api/python/)
- [AWS Glue Pricing](https://aws.amazon.com/glue/pricing/)
