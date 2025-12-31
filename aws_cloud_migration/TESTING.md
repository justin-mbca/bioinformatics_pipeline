# Testing Guide

## Quick Start

### Local Testing (No AWS Account)
```bash
cd aws_cloud_migration
bash scripts/run_local_tests.sh
```

### AWS Integration Testing
```bash
# Configure AWS CLI first
aws configure

# Run AWS tests
bash scripts/run_aws_tests.sh
```

## Test Coverage
- Healthcare data standards (FHIR, HL7, OMOP)
- Docker containers
- Terraform configurations
- Python code quality
- AWS services integration

## CI/CD
Tests run automatically on every push via GitHub Actions.
