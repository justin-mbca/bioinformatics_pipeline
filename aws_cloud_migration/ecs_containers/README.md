# ECS Containers for RNA-Seq Analysis

## Overview

This directory contains Docker containers and ECS/Fargate configurations for running RNA-Seq analysis workflows.

## Files

### Dockerfile.rnaseq
Production-ready Docker image with:
- Python 3.9
- Bioinformatics packages (numpy, pandas, scipy, biopython, HTSeq, pysam)
- AWS SDK (boto3)
- Non-root user for security
- Health checks

### ecs_task_definition.json
Fargate task definition with:
- 2 vCPU, 8 GB memory
- CloudWatch logging
- Environment variables
- IAM roles

### fargate_deployment.yaml
CloudFormation template for:
- ECS cluster
- ECS service
- Auto-scaling policies
- Security groups
- IAM roles

## Building the Image

```bash
# Build locally
docker build -t bioinformatics-rnaseq -f Dockerfile.rnaseq .

# Test locally
docker run -it bioinformatics-rnaseq python -c "import numpy; print(numpy.__version__)"

# Push to ECR
aws ecr create-repository --repository-name bioinformatics-rnaseq
aws ecr get-login-password | docker login --username AWS --password-stdin <account>.dkr.ecr.<region>.amazonaws.com
docker tag bioinformatics-rnaseq:latest <account>.dkr.ecr.<region>.amazonaws.com/bioinformatics-rnaseq:latest
docker push <account>.dkr.ecr.<region>.amazonaws.com/bioinformatics-rnaseq:latest
```

## Deployment

### Option 1: AWS CLI
```bash
# Register task definition
aws ecs register-task-definition --cli-input-json file://ecs_task_definition.json

# Create cluster
aws ecs create-cluster --cluster-name bioinformatics-cluster

# Run task
aws ecs run-task \
    --cluster bioinformatics-cluster \
    --task-definition bioinformatics-rnaseq-analysis \
    --launch-type FARGATE
```

### Option 2: CloudFormation
```bash
aws cloudformation create-stack \
    --stack-name bioinformatics-ecs \
    --template-body file://fargate_deployment.yaml \
    --capabilities CAPABILITY_IAM
```

## Cost

- Fargate: ~$0.04/vCPU-hour + $0.004/GB-hour
- Example: 2 vCPU, 8 GB for 1 hour = $0.11

## Security

- Non-root user in container
- Read-only root filesystem
- Security groups restrict network access
- IAM roles for AWS service access
- Secrets via AWS Secrets Manager

## Additional Resources

- [Amazon ECS Documentation](https://docs.aws.amazon.com/ecs/)
- [Docker Best Practices](https://docs.docker.com/develop/dev-best-practices/)
