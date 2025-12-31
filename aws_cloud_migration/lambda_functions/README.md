# AWS Lambda Functions

## Overview

This directory contains AWS Lambda functions for pipeline automation and data quality validation.

## Functions

### 1. trigger_pipeline.py
**Purpose**: Automatically triggers AWS Glue ETL job when new RNA-Seq files are uploaded to S3.

**Features**:
- S3 event-driven execution
- Validates file format and size
- Starts Glue job with parameters
- Returns job run ID
- Error handling and logging

**Environment Variables**:
- `GLUE_JOB_NAME`: Name of Glue job to trigger
- `GLUE_DATABASE`: Glue database name
- `GLUE_TABLE`: Glue table name
- `OUTPUT_BUCKET`: S3 bucket for output
- `COUNT_THRESHOLD`: Minimum count threshold

### 2. data_quality_check.py
**Purpose**: Validates RNA-Seq data files for quality issues and sends notifications.

**Features**:
- File format validation (CSV/TSV)
- Column header validation
- Missing value detection
- Numeric data validation
- SNS notifications on failure

**Environment Variables**:
- `SNS_TOPIC_ARN`: SNS topic for notifications
- `MAX_MISSING_PERCENT`: Maximum allowed missing values (%)
- `MIN_SAMPLE_COUNT`: Minimum number of samples required
- `REQUIRED_COLUMNS`: Comma-separated list of required columns

## Deployment

### Option 1: SAM CLI
```bash
sam build
sam deploy --guided
```

### Option 2: AWS CLI
```bash
# Package Lambda function
zip -r trigger_pipeline.zip trigger_pipeline.py
aws lambda create-function \
    --function-name trigger-pipeline \
    --runtime python3.9 \
    --role arn:aws:iam::ACCOUNT:role/lambda-role \
    --handler trigger_pipeline.lambda_handler \
    --zip-file fileb://trigger_pipeline.zip
```

### Option 3: Terraform
See `../terraform/lambda.tf`

## Testing

### Unit Tests
```bash
python -m pytest tests/
```

### Local Testing
```bash
# Test with sample event
python trigger_pipeline.py
python data_quality_check.py
```

## Monitoring

- CloudWatch Logs: `/aws/lambda/<function-name>`
- CloudWatch Metrics: Invocations, Errors, Duration
- X-Ray: Distributed tracing

## Cost

- $0.20 per 1M requests
- $0.0000166667 per GB-second
- First 1M requests/month free

## Additional Resources

- [AWS Lambda Documentation](https://docs.aws.amazon.com/lambda/)
- [Lambda Best Practices](https://docs.aws.amazon.com/lambda/latest/dg/best-practices.html)
