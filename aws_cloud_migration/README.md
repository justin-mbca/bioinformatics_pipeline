# AWS Cloud Migration - Bioinformatics Pipeline

## Overview

This directory contains production-ready AWS cloud infrastructure and healthcare data standards integration for the bioinformatics pipeline. The implementation demonstrates expertise with AWS services (Glue, EMR, Lambda, S3, ECS) and healthcare data formats (FHIR, HL7, OMOP) suitable for Medical Data Specialist roles.

## Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                        AWS Cloud Infrastructure                  │
└─────────────────────────────────────────────────────────────────┘

 ┌──────────────┐
 │ S3 Data Lake │ ◄────────────┐
 └──────┬───────┘              │
        │                      │
        ▼                      │
 ┌──────────────┐       ┌─────┴──────┐      ┌────────────────┐
 │  AWS Glue    │       │  Lambda    │      │  Healthcare    │
 │  ETL Jobs    │◄─────▶│  Functions │◄────▶│  Standards     │
 └──────┬───────┘       └────────────┘      │  (FHIR/HL7)    │
        │                                    └────────────────┘
        ▼
 ┌──────────────┐                           ┌────────────────┐
 │  AWS EMR     │                           │  ECS/Fargate   │
 │  Spark Jobs  │                           │  Containers    │
 └──────────────┘                           └────────────────┘
```

## Directory Structure

```
aws_cloud_migration/
├── README.md                          # This file
├── glue_jobs/                         # AWS Glue ETL jobs
│   ├── etl_rnaseq_to_s3.py           # Main ETL script
│   ├── glue_job_definition.json      # Job configuration
│   ├── data_catalog_setup.py         # Catalog creation
│   └── README.md                      # Glue documentation
├── emr_spark/                         # EMR Spark analysis
│   ├── spark_deseq2_analysis.py      # DESeq2-style analysis
│   ├── emr_cluster_config.json       # Cluster configuration
│   ├── bootstrap_actions.sh          # Setup script
│   └── README.md                      # EMR documentation
├── lambda_functions/                  # AWS Lambda functions
│   ├── trigger_pipeline.py           # Pipeline automation
│   ├── data_quality_check.py         # Data validation
│   ├── lambda_deployment.yaml        # SAM template
│   └── README.md                      # Lambda documentation
├── s3_data_lake/                      # S3 bucket structure
│   ├── lifecycle_policies.json       # Storage policies
│   ├── upload_datasets.py            # Upload script
│   └── README.md                      # S3 documentation
├── ecs_containers/                    # ECS containers
│   ├── Dockerfile.rnaseq             # Container image
│   ├── ecs_task_definition.json      # Task config
│   ├── fargate_deployment.yaml       # Service config
│   └── README.md                      # ECS documentation
├── terraform/                         # Infrastructure as Code
│   ├── main.tf                       # Provider config
│   ├── variables.tf                  # Input variables
│   ├── s3.tf                        # S3 resources
│   ├── glue.tf                      # Glue resources
│   ├── emr.tf                       # EMR resources
│   ├── lambda.tf                    # Lambda resources
│   ├── outputs.tf                   # Output values
│   └── README.md                     # Terraform docs
└── healthcare_standards/             # Healthcare integration
    ├── fhir_genomics_integration.py  # FHIR R4 implementation
    ├── hl7_v2_parser.py             # HL7 message parser
    ├── omop_cdm_mapping.py          # OMOP CDM mapping
    └── README.md                     # Healthcare docs
```

## AWS Services Used

### 1. AWS Glue
- **Purpose**: Serverless ETL for RNA-Seq data processing
- **Features**: Data catalog, schema discovery, Spark-based transformations
- **Cost**: ~$0.44/DPU-hour
- **Documentation**: [glue_jobs/README.md](glue_jobs/README.md)

### 2. AWS EMR (Elastic MapReduce)
- **Purpose**: Large-scale differential expression analysis with Spark
- **Features**: Auto-scaling, managed Hadoop/Spark, spot instance support
- **Cost**: ~$0.96/hour (m5.xlarge + 3x r5.xlarge)
- **Documentation**: [emr_spark/README.md](emr_spark/README.md)

### 3. AWS Lambda
- **Purpose**: Event-driven pipeline automation and data validation
- **Features**: S3 triggers, SNS notifications, automatic scaling
- **Cost**: ~$0.20 per 1M requests
- **Documentation**: [lambda_functions/README.md](lambda_functions/README.md)

### 4. Amazon S3
- **Purpose**: Scalable data lake for genomics datasets
- **Features**: Lifecycle policies, versioning, encryption, analytics
- **Cost**: ~$0.023/GB/month (Standard), less for IA/Glacier
- **Documentation**: [s3_data_lake/README.md](s3_data_lake/README.md)

### 5. AWS ECS/Fargate
- **Purpose**: Containerized analysis workflows
- **Features**: Serverless containers, auto-scaling, CloudWatch integration
- **Cost**: ~$0.04/vCPU-hour + $0.004/GB-hour
- **Documentation**: [ecs_containers/README.md](ecs_containers/README.md)

## Healthcare Data Standards

### FHIR R4 (Fast Healthcare Interoperability Resources)
- Create genomic observations with LOINC codes
- Support for gene expression values (TPM, FPKM, counts)
- Patient and specimen references
- JSON serialization for EHR integration

### HL7 v2.x Message Parsing
- Parse ORU^R01 (Lab Results) messages
- Extract patient demographics and observations
- Map to RNA-Seq metadata format
- Error handling for malformed messages

### OMOP Common Data Model
- Map patient demographics to PERSON table
- Map tissue samples to SPECIMEN table
- Map gene expression to MEASUREMENT table
- Generate SQL INSERT statements

**Documentation**: [healthcare_standards/README.md](healthcare_standards/README.md)

## Quick Start

### Prerequisites

1. **AWS Account**: Active AWS account with appropriate permissions
2. **AWS CLI**: Configured with credentials
   ```bash
   aws configure
   ```
3. **Terraform**: For infrastructure deployment
   ```bash
   terraform version  # Should be >= 1.0
   ```

### Option 1: Deploy with Terraform (Recommended)

```bash
# Navigate to terraform directory
cd terraform/

# Initialize Terraform
terraform init

# Review the deployment plan
terraform plan

# Deploy infrastructure
terraform apply

# View outputs
terraform output
```

### Option 2: Manual Deployment

Follow the detailed setup instructions in each component's README:

1. **S3 Data Lake**: [s3_data_lake/README.md](s3_data_lake/README.md)
2. **AWS Glue**: [glue_jobs/README.md](glue_jobs/README.md)
3. **AWS EMR**: [emr_spark/README.md](emr_spark/README.md)
4. **AWS Lambda**: [lambda_functions/README.md](lambda_functions/README.md)
5. **ECS Containers**: [ecs_containers/README.md](ecs_containers/README.md)

## Workflow Examples

### Example 1: Automated RNA-Seq Processing

```
1. Upload raw counts → S3 (raw/counts/)
2. Lambda triggers Glue ETL job
3. Glue processes and normalizes data
4. Results stored in S3 (processed/normalized/)
5. SNS notification sent on completion
```

### Example 2: Large-Scale Differential Expression

```
1. Submit Spark job to EMR
2. Read normalized data from S3
3. Calculate DESeq2-style statistics
4. Write results to S3 (results/de_genes/)
5. Export significant genes to FHIR format
```

### Example 3: Healthcare System Integration

```
1. Receive HL7 message from EHR
2. Parse patient demographics and lab results
3. Map to RNA-Seq metadata format
4. Trigger analysis pipeline
5. Export results to OMOP CDM database
```

## Cost Optimization

### Monthly Cost Estimates (us-east-1)

| Service | Usage | Monthly Cost |
|---------|-------|--------------|
| S3 (1 TB data) | Standard + IA + Glacier | $23 |
| Glue (10 DPU-hrs/day) | ETL processing | $132 |
| EMR (8 hrs/day) | Spark analysis | $230 |
| Lambda (1M requests) | Pipeline automation | $5 |
| ECS Fargate (100 hrs) | Container tasks | $20 |
| Data Transfer | Within AWS | $50 |
| **Total** | | **~$460/month** |

### Optimization Strategies

1. **Use Spot Instances**: 50-70% savings on EMR
2. **S3 Lifecycle Policies**: Auto-transition to cheaper storage
3. **EMR Auto-Termination**: Prevent idle cluster charges
4. **Lambda Reserved Concurrency**: Predictable costs
5. **S3 Intelligent-Tiering**: Automatic cost optimization
6. **Glue Job Bookmarks**: Process only new data

## Security Best Practices

### Data Protection
- **Encryption at Rest**: AES-256 for S3, KMS for EBS
- **Encryption in Transit**: TLS 1.2+ for all connections
- **Bucket Policies**: Deny insecure transport
- **Versioning**: Protect against accidental deletion

### Access Control
- **IAM Roles**: Least privilege access
- **Resource-Based Policies**: Fine-grained permissions
- **VPC Endpoints**: Private S3 access
- **Security Groups**: Restrict network access

### Compliance
- **HIPAA**: Enable logging, encryption, access controls
- **GDPR**: Data subject rights, data minimization
- **FDA 21 CFR Part 11**: Electronic signatures, audit trails
- **SOC 2**: Security and availability controls

### Monitoring
- **CloudTrail**: API audit logging
- **CloudWatch**: Metrics and alarms
- **VPC Flow Logs**: Network traffic analysis
- **AWS Config**: Resource compliance tracking

## Performance Tuning

### AWS Glue
- Use G.2X workers for memory-intensive jobs
- Enable job bookmarks for incremental processing
- Partition large datasets by date/chromosome
- Use Parquet format with Snappy compression

### AWS EMR
- Use r5 instances for memory-intensive workloads
- Enable dynamic allocation and adaptive execution
- Cache frequently accessed data
- Use spot instances for task nodes

### AWS Lambda
- Increase memory for CPU-intensive tasks
- Use environment variables for configuration
- Implement connection pooling for databases
- Enable X-Ray for distributed tracing

### Amazon S3
- Use multipart uploads for large files
- Enable transfer acceleration for global uploads
- Use S3 Select for filtering data
- Implement request rate optimization

## Troubleshooting

### Common Issues

**Glue Job Failures**
- Check CloudWatch logs: `/aws-glue/jobs/output`
- Verify IAM role permissions
- Ensure S3 paths are correct
- Increase DPU allocation if OOM

**EMR Cluster Won't Start**
- Verify subnet has internet access
- Check EC2 service limits
- Ensure IAM roles exist
- Review bootstrap action logs

**Lambda Timeouts**
- Increase function timeout (max 15 min)
- Optimize code performance
- Use asynchronous invocations
- Check network connectivity

**S3 Access Denied**
- Verify bucket policy
- Check IAM role permissions
- Ensure bucket exists in correct region
- Review VPC endpoint configuration

## Monitoring and Alerts

### Key Metrics to Monitor

**AWS Glue**
- Job success/failure rate
- DPU hours consumed
- Data processed (GB)
- Job execution time

**AWS EMR**
- Cluster utilization (CPU, memory)
- HDFS/YARN metrics
- Step success rate
- Node failures

**AWS Lambda**
- Invocation count
- Error rate
- Duration
- Throttles

**Amazon S3**
- Request rate
- Error rate (4xx, 5xx)
- Bytes downloaded/uploaded
- Bucket size

### CloudWatch Alarms

Set up alarms for:
- Glue job failures
- EMR cluster errors
- Lambda error rate > 1%
- S3 bucket size > threshold
- Unusual access patterns

## Best Practices

### Development Workflow
1. Test locally with sample data
2. Deploy to dev environment
3. Run integration tests
4. Deploy to staging
5. Perform UAT
6. Deploy to production

### Data Management
1. Use consistent naming conventions
2. Implement data versioning
3. Document data lineage
4. Regular data quality checks
5. Automated backup and recovery

### Code Quality
1. Follow PEP 8 style guide
2. Add type hints to functions
3. Write comprehensive docstrings
4. Include error handling
5. Add logging statements

### Infrastructure
1. Use Terraform for all resources
2. Version control configurations
3. Implement CI/CD pipelines
4. Regular security audits
5. Cost monitoring and optimization

## Additional Resources

### AWS Documentation
- [AWS Glue Developer Guide](https://docs.aws.amazon.com/glue/)
- [Amazon EMR Management Guide](https://docs.aws.amazon.com/emr/)
- [AWS Lambda Developer Guide](https://docs.aws.amazon.com/lambda/)
- [Amazon S3 User Guide](https://docs.aws.amazon.com/s3/)
- [Amazon ECS Developer Guide](https://docs.aws.amazon.com/ecs/)

### Healthcare Standards
- [FHIR Genomics](http://hl7.org/fhir/uv/genomics-reporting/)
- [HL7 v2 Specification](http://www.hl7.org/implement/standards/)
- [OMOP CDM Documentation](https://ohdsi.github.io/CommonDataModel/)

### Best Practices
- [AWS Well-Architected Framework](https://aws.amazon.com/architecture/well-architected/)
- [AWS Security Best Practices](https://aws.amazon.com/security/best-practices/)
- [Terraform Best Practices](https://www.terraform-best-practices.com/)

## Support and Contributing

For questions, issues, or contributions:
- Open an issue in the repository
- Contact the bioinformatics team
- Review the component-specific README files

## License

This project is licensed under the MIT License. See the LICENSE file for details.

## Acknowledgments

This implementation demonstrates production-ready patterns for:
- Cloud-native bioinformatics workflows
- Healthcare data standards integration
- Enterprise security and compliance
- Cost-optimized infrastructure design
