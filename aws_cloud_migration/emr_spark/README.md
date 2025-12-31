# AWS EMR Spark Jobs

This directory contains Apache Spark jobs for large-scale differential expression analysis on AWS EMR (Elastic MapReduce).

## Overview

AWS EMR is a cloud-native big data platform that makes it easy to process vast amounts of data using Apache Spark, Hadoop, and other distributed frameworks. These scripts implement DESeq2-style differential expression analysis at scale.

## Files

### `spark_deseq2_analysis.py`
Production-grade PySpark implementation of differential expression analysis featuring:
- **Size Factor Calculation**: DESeq2-style median-of-ratios normalization
- **Log2 Fold Change**: Calculation between treatment and control groups
- **Statistical Testing**: Wald test approximation with FDR correction
- **Scalability**: Processes datasets with millions of genes and thousands of samples
- **Fault Tolerance**: Built on Spark's resilient distributed datasets (RDDs)

**Key Features:**
- Type-hinted class-based design
- Comprehensive logging and error handling
- Optimized for S3 data I/O with Parquet format
- Adaptive query execution for performance
- Supports dynamic resource allocation

### `emr_cluster_config.json`
Complete EMR cluster configuration including:
- **Release**: EMR 6.10.0 with Spark 3.3.1
- **Applications**: Spark, Hadoop, Hive, Livy, JupyterHub
- **Master Node**: m5.xlarge (4 vCPU, 16 GB RAM)
- **Core Nodes**: r5.xlarge (4 vCPU, 32 GB RAM) - memory optimized for genomics
- **Auto-scaling**: 2-10 instances based on YARN memory utilization
- **Storage**: GP3 EBS volumes with 3000 IOPS
- **Auto-termination**: Cluster terminates after 1 hour of inactivity

**Cost Optimization:**
- Managed scaling policy for dynamic resource allocation
- Auto-termination to prevent idle charges
- Step concurrency for parallel job execution
- Spot instances option (modify "Market" to "SPOT")

### `bootstrap_actions.sh`
Comprehensive bootstrap script that:
- Installs bioinformatics Python packages (numpy, pandas, scipy, biopython, HTSeq, pysam)
- Configures Spark settings for genomics workloads
- Sets up logging and monitoring
- Optimizes memory settings for large datasets
- Downloads custom analysis scripts from S3

**Optimizations:**
- Increased network timeouts for large S3 transfers
- Kryo serialization for better performance
- Adaptive query execution enabled
- Memory overhead configured for large shuffles

## Prerequisites

### 1. AWS Account Setup
- Active AWS account with EMR permissions
- AWS CLI configured with credentials

### 2. IAM Roles
Create two IAM roles:

**EMR Service Role** (`EMR_DefaultRole`):
```json
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Effect": "Allow",
      "Principal": {
        "Service": "elasticmapreduce.amazonaws.com"
      },
      "Action": "sts:AssumeRole"
    }
  ]
}
```

**EMR EC2 Instance Profile** (`EMR_EC2_DefaultRole`):
- S3 read/write access
- CloudWatch Logs write access
- EC2 permissions for cluster communication

### 3. S3 Bucket Structure
```
s3://bioinformatics-pipeline/
├── processed/normalized/         # Input: Processed count data
├── results/de_genes/             # Output: DE analysis results
├── logs/emr/                     # EMR cluster logs
├── bootstrap/
│   └── bootstrap_actions.sh      # Bootstrap script
└── emr-scripts/
    └── spark_deseq2_analysis.py  # Analysis script
```

### 4. EC2 Key Pair
Create an EC2 key pair for SSH access:
```bash
aws ec2 create-key-pair --key-name emr-key-pair --query 'KeyMaterial' --output text > emr-key-pair.pem
chmod 400 emr-key-pair.pem
```

### 5. VPC and Subnet
- VPC with internet access (or NAT gateway)
- Public or private subnet for EMR cluster
- Security groups allowing EMR communication

## Setup Instructions

### 1. Upload Scripts to S3

```bash
# Upload bootstrap script
aws s3 cp bootstrap_actions.sh s3://bioinformatics-pipeline/bootstrap/

# Upload analysis script
aws s3 cp spark_deseq2_analysis.py s3://bioinformatics-pipeline/emr-scripts/
```

### 2. Update Configuration

Edit `emr_cluster_config.json`:
- Replace `subnet-xxxxxxxx` with your subnet ID
- Update `Ec2KeyName` with your key pair name
- Adjust instance types and counts as needed

### 3. Create EMR Cluster

Using AWS CLI:
```bash
aws emr create-cluster --cli-input-json file://emr_cluster_config.json
```

Using Terraform (see `../terraform/emr.tf`):
```bash
cd ../terraform
terraform init
terraform plan
terraform apply
```

### 4. Submit Spark Job

**Option A: Add Step to Running Cluster**

```bash
CLUSTER_ID="j-XXXXXXXXXXXXX"

aws emr add-steps \
    --cluster-id $CLUSTER_ID \
    --steps Type=Spark,Name="DESeq2 Analysis",ActionOnFailure=CONTINUE,Args=[
        --deploy-mode,cluster,
        --master,yarn,
        --conf,spark.dynamicAllocation.enabled=true,
        --conf,spark.sql.adaptive.enabled=true,
        s3://bioinformatics-pipeline/emr-scripts/spark_deseq2_analysis.py,
        s3://bioinformatics-pipeline/processed/normalized/,
        s3://bioinformatics-pipeline/results/de_genes/,
        'sample_01,sample_02,sample_03',
        'sample_04,sample_05,sample_06'
    ]
```

**Option B: Submit via spark-submit (SSH to master node)**

```bash
# SSH to master node
ssh -i emr-key-pair.pem hadoop@<master-public-dns>

# Submit job
spark-submit \
    --master yarn \
    --deploy-mode cluster \
    --num-executors 6 \
    --executor-memory 8G \
    --executor-cores 2 \
    --driver-memory 4G \
    --conf spark.dynamicAllocation.enabled=true \
    --conf spark.sql.adaptive.enabled=true \
    s3://bioinformatics-pipeline/emr-scripts/spark_deseq2_analysis.py \
    s3://bioinformatics-pipeline/processed/normalized/ \
    s3://bioinformatics-pipeline/results/de_genes/ \
    'sample_01,sample_02,sample_03' \
    'sample_04,sample_05,sample_06'
```

## Usage Examples

### Basic Usage

```bash
spark-submit spark_deseq2_analysis.py \
    <input_path> \
    <output_path> \
    <condition1_samples> \
    <condition2_samples>
```

### Example: Treatment vs Control

```bash
spark-submit spark_deseq2_analysis.py \
    s3://bioinformatics-pipeline/processed/normalized/ \
    s3://bioinformatics-pipeline/results/de_genes/treatment_vs_control/ \
    'treated_01,treated_02,treated_03,treated_04' \
    'control_01,control_02,control_03,control_04'
```

### Example: Time Point Comparison

```bash
spark-submit spark_deseq2_analysis.py \
    s3://bioinformatics-pipeline/processed/normalized/ \
    s3://bioinformatics-pipeline/results/de_genes/day7_vs_day0/ \
    'day7_rep1,day7_rep2,day7_rep3' \
    'day0_rep1,day0_rep2,day0_rep3'
```

## Monitoring

### 1. EMR Console
Monitor cluster status:
1. AWS Console → EMR → Clusters
2. Select your cluster
3. View steps, hardware, and application metrics

### 2. Spark UI
Access Spark UI for job details:
1. Enable SSH tunnel to master node:
   ```bash
   aws emr socks --cluster-id $CLUSTER_ID --key-pair-file emr-key-pair.pem
   ```
2. Configure browser proxy (localhost:8157)
3. Navigate to master node public DNS:18080

### 3. CloudWatch Logs
View logs in CloudWatch:
```bash
aws logs tail /aws/emr/$CLUSTER_ID --follow
```

### 4. Application Logs
View Spark application logs:
```bash
# On master node
yarn logs -applicationId <application_id>
```

## Output Format

The analysis produces two output directories:

### `all_genes/`
Parquet files containing results for all genes:
```
gene_id | gene_name | chromosome | baseMean | log2FoldChange | lfcSE | pvalue | padj | cond1_mean | cond2_mean
```

### `significant_genes/`
Filtered results (padj < 0.05, |log2FC| > 1.0):
- Contains only significantly differentially expressed genes
- Sorted by adjusted p-value

**Query Results with Athena:**
```sql
CREATE EXTERNAL TABLE de_results (
    gene_id STRING,
    gene_name STRING,
    baseMean DOUBLE,
    log2FoldChange DOUBLE,
    pvalue DOUBLE,
    padj DOUBLE
)
STORED AS PARQUET
LOCATION 's3://bioinformatics-pipeline/results/de_genes/treatment_vs_control/all_genes/';

SELECT gene_name, log2FoldChange, padj
FROM de_results
WHERE padj < 0.01
ORDER BY ABS(log2FoldChange) DESC
LIMIT 20;
```

## Performance Tuning

### For Small Datasets (< 10 GB)
```json
{
  "InstanceType": "m5.xlarge",
  "InstanceCount": 2,
  "MaxCapacity": 4
}
```

### For Medium Datasets (10-100 GB)
```json
{
  "InstanceType": "r5.xlarge",
  "InstanceCount": 3,
  "MaxCapacity": 10
}
```

### For Large Datasets (> 100 GB)
```json
{
  "InstanceType": "r5.2xlarge",
  "InstanceCount": 5,
  "MaxCapacity": 20
}
```

**Spark Configuration for Large Data:**
```bash
--executor-memory 16G \
--executor-cores 4 \
--driver-memory 8G \
--conf spark.sql.shuffle.partitions=400 \
--conf spark.default.parallelism=400
```

## Cost Estimates

**Example Configuration** (us-east-1):
- Master: 1 x m5.xlarge = $0.192/hr
- Core: 3 x r5.xlarge = $0.252/hr each = $0.756/hr
- Total: $0.948/hr

**Processing Time Estimates:**
- 100 genes, 10 samples: ~5 minutes = $0.08
- 20,000 genes, 100 samples: ~30 minutes = $0.47
- 50,000 genes, 1000 samples: ~2 hours = $1.90

**Cost Optimization:**
- Use Spot instances for 50-70% savings
- Enable auto-termination
- Use managed scaling
- Process multiple datasets in parallel (step concurrency)

## Troubleshooting

### Job Fails with OutOfMemoryError

Increase executor memory:
```bash
--executor-memory 16G \
--conf spark.executor.memoryOverhead=4096
```

### Slow S3 Read/Write

Increase parallelism:
```bash
--conf spark.sql.files.maxPartitionBytes=128m \
--conf spark.sql.shuffle.partitions=400
```

### Cluster Won't Start

Check:
- IAM roles have correct permissions
- Subnet has internet access or NAT gateway
- Security groups allow EMR communication
- EC2 limits not exceeded

### Bootstrap Action Fails

View bootstrap logs:
```bash
aws s3 ls s3://bioinformatics-pipeline/logs/emr/$CLUSTER_ID/node/
aws s3 cp s3://bioinformatics-pipeline/logs/emr/$CLUSTER_ID/node/$NODE_ID/bootstrap-actions/ . --recursive
```

## Best Practices

1. **Data Partitioning**: Partition large datasets by chromosome or experimental batch
2. **Checkpointing**: Enable Spark checkpointing for long-running jobs
3. **Caching**: Cache DataFrames that are reused multiple times
4. **Resource Allocation**: Match instance types to workload (CPU vs memory)
5. **Logging**: Enable detailed logging for production jobs
6. **Testing**: Test with small datasets before running on full data
7. **Security**: Use VPC, encryption, and IAM roles properly

## Security Considerations

- **Encryption**: Enable S3 encryption at rest and in transit
- **Network**: Run EMR in private subnet with NAT gateway
- **Access**: Use IAM roles instead of access keys
- **Logging**: Enable CloudTrail for audit logs
- **Security Groups**: Restrict access to master node
- **Data**: Comply with HIPAA/GDPR requirements for sensitive data

## Integration

### With AWS Glue
Read from Glue Data Catalog:
```python
df = spark.read.format("parquet").table("rnaseq_db.processed_counts")
```

### With Amazon Athena
Query results directly from S3 without loading into memory

### With SageMaker
Export significant genes for machine learning model training

## Additional Resources

- [AWS EMR Documentation](https://docs.aws.amazon.com/emr/)
- [Apache Spark Documentation](https://spark.apache.org/docs/latest/)
- [EMR Best Practices](https://docs.aws.amazon.com/emr/latest/ManagementGuide/emr-plan-instances-guidelines.html)
- [Spark Performance Tuning](https://spark.apache.org/docs/latest/sql-performance-tuning.html)
- [EMR Pricing](https://aws.amazon.com/emr/pricing/)
