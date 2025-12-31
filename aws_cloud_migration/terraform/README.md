# Terraform Infrastructure as Code

## Overview

This directory contains Terraform configurations to deploy the complete AWS infrastructure for the bioinformatics pipeline.

## Prerequisites

1. **Terraform**: Install Terraform >= 1.0
   ```bash
   curl -fsSL https://apt.releases.hashicorp.com/gpg | sudo apt-key add -
   sudo apt-add-repository "deb [arch=amd64] https://apt.releases.hashicorp.com $(lsb_release -cs) main"
   sudo apt-get update && sudo apt-get install terraform
   ```

2. **AWS CLI**: Configure AWS credentials
   ```bash
   aws configure
   ```

3. **S3 Backend**: Create S3 bucket and DynamoDB table for state
   ```bash
   aws s3 mb s3://bioinformatics-terraform-state
   aws dynamodb create-table \
       --table-name terraform-state-lock \
       --attribute-definitions AttributeName=LockID,AttributeType=S \
       --key-schema AttributeName=LockID,KeyType=HASH \
       --billing-mode PAY_PER_REQUEST
   ```

## Quick Start

1. **Initialize Terraform**:
   ```bash
   terraform init
   ```

2. **Review the plan**:
   ```bash
   terraform plan
   ```

3. **Apply the configuration**:
   ```bash
   terraform apply
   ```

4. **Destroy resources** (when needed):
   ```bash
   terraform destroy
   ```

## Configuration Files

| File | Purpose |
|------|---------|
| `main.tf` | Provider and backend configuration |
| `variables.tf` | Input variables |
| `outputs.tf` | Output values |
| `s3.tf` | S3 data lake bucket |
| `glue.tf` | AWS Glue resources |
| `emr.tf` | EMR cluster |
| `lambda.tf` | Lambda functions |

## Customization

Create a `terraform.tfvars` file to customize variables:

```hcl
aws_region          = "us-east-1"
environment         = "production"
s3_bucket_name      = "my-bioinformatics-pipeline"
emr_instance_type   = "r5.xlarge"
sns_email           = "my-email@example.com"
```

## Resources Created

- **S3 Bucket**: Data lake with lifecycle policies
- **Glue Database**: Data catalog
- **Glue ETL Job**: RNA-Seq processing
- **Glue Crawler**: Automatic schema discovery
- **EMR Cluster**: Spark analysis cluster
- **Lambda Functions**: Pipeline automation
- **SNS Topic**: Notifications
- **IAM Roles**: Service permissions
- **Security Groups**: Network access control
- **CloudWatch Logs**: Monitoring and logging

## Cost Estimates (us-east-1)

**Monthly Costs** (approximate):
- S3 (1 TB): ~$23
- Glue (10 DPU-hours/day): ~$132
- EMR (running 8 hrs/day): ~$230
- Lambda (1M invocations): ~$5
- Data Transfer: ~$50
- **Total**: ~$440/month

**Optimization Tips**:
- Use Spot instances for EMR (50-70% savings)
- Enable EMR auto-termination
- Use S3 lifecycle policies
- Schedule Glue crawlers appropriately

## Security Best Practices

1. **Encryption**: All data encrypted at rest and in transit
2. **IAM**: Least privilege access policies
3. **VPC**: Resources in private subnets
4. **Logging**: CloudTrail and CloudWatch enabled
5. **Secrets**: Use AWS Secrets Manager for sensitive data

## Outputs

After successful `terraform apply`, you'll see:

```
s3_bucket_name         = "bioinformatics-pipeline"
glue_database_name     = "rnaseq_db"
emr_cluster_id         = "j-XXXXXXXXXXXXX"
trigger_pipeline_lambda_arn = "arn:aws:lambda:..."
```

## Troubleshooting

**State Lock Issues**:
```bash
terraform force-unlock <LOCK_ID>
```

**Resource Already Exists**:
```bash
terraform import <resource_type>.<name> <resource_id>
```

**Permission Denied**:
- Check IAM permissions
- Verify AWS CLI configuration

## Additional Resources

- [Terraform AWS Provider Docs](https://registry.terraform.io/providers/hashicorp/aws/latest/docs)
- [AWS Terraform Best Practices](https://aws.amazon.com/blogs/apn/terraform-beyond-the-basics-with-aws/)
- [Terraform State Management](https://www.terraform.io/docs/language/state/index.html)
