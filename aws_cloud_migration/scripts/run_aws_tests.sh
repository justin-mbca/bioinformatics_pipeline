#!/bin/bash
# AWS integration tests - requires AWS credentials
# Tests: S3, Lambda, Glue, EMR, ECS

set -e

echo "☁️  Starting AWS Integration Tests..."

# Check AWS credentials
if ! aws sts get-caller-identity &>/dev/null; then
    echo "❌ AWS credentials not configured. Run 'aws configure' first."
    exit 1
fi

# Get AWS account ID
AWS_ACCOUNT_ID=$(aws sts get-caller-identity --query Account --output text)
AWS_REGION=${AWS_REGION:-us-east-1}

echo "Using AWS Account: $AWS_ACCOUNT_ID"
echo "Using AWS Region: $AWS_REGION"

# Generate unique names
TIMESTAMP=$(date +%s)
TEST_BUCKET="bioinformatics-test-${AWS_ACCOUNT_ID}-${TIMESTAMP}"
TEST_PREFIX="test-${TIMESTAMP}"

# Cleanup function
cleanup() {
    echo "🧹 Cleaning up AWS resources..."
    
    # Delete S3 bucket
    aws s3 rb s3://${TEST_BUCKET} --force 2>/dev/null || true
    
    # Delete Lambda function
    aws lambda delete-function --function-name ${TEST_PREFIX}-trigger 2>/dev/null || true
    
    # Delete Glue database
    aws glue delete-database --name ${TEST_PREFIX}_db 2>/dev/null || true
    
    echo "✓ Cleanup complete"
}

trap cleanup EXIT

# Test 1: S3 Operations
echo -e "\n🪣 Test 1: S3 Data Lake"
aws s3 mb s3://${TEST_BUCKET}
echo "test data" > test_file.txt
aws s3 cp test_file.txt s3://${TEST_BUCKET}/raw/test_file.txt
aws s3 ls s3://${TEST_BUCKET}/raw/
echo "✓ S3 test passed"

# Test 2: S3 Upload Script
echo -e "\n📤 Test 2: S3 Upload Script"
cd s3_data_lake
python3 upload_datasets.py \
    --bucket ${TEST_BUCKET} \
    --file ../test_file.txt \
    --key test/upload_script_test.txt
cd ..
echo "✓ Upload script test passed"

# Test 3: Lambda Function Deployment
echo -e "\n⚡ Test 3: Lambda Function"
cd lambda_functions
zip -q function.zip trigger_pipeline.py

# Create IAM role (if doesn't exist)
ROLE_NAME="${TEST_PREFIX}-lambda-role"
aws iam create-role --role-name ${ROLE_NAME} \
    --assume-role-policy-document '{
        "Version": "2012-10-17",
        "Statement": [{
            "Effect": "Allow",
            "Principal": {"Service": "lambda.amazonaws.com"},
            "Action": "sts:AssumeRole"
        }]
    }' 2>/dev/null || true

aws iam attach-role-policy --role-name ${ROLE_NAME} \
    --policy-arn arn:aws:iam::aws:policy/service-role/AWSLambdaBasicExecutionRole 2>/dev/null || true

# Wait for role propagation
sleep 10

ROLE_ARN="arn:aws:iam::${AWS_ACCOUNT_ID}:role/${ROLE_NAME}"

aws lambda create-function \
    --function-name ${TEST_PREFIX}-trigger \
    --runtime python3.9 \
    --role ${ROLE_ARN} \
    --handler trigger_pipeline.lambda_handler \
    --zip-file fileb://function.zip \
    --timeout 60

# Test invocation
aws lambda invoke \
    --function-name ${TEST_PREFIX}-trigger \
    --payload '{"Records":[{"s3":{"bucket":{"name":"'${TEST_BUCKET}'"},"object":{"key":"test.csv"}}}]}' \
    response.json

cat response.json
cd ..
echo "✓ Lambda test passed"

# Test 4: Glue Data Catalog
echo -e "\n📊 Test 4: Glue Data Catalog"
cd glue_jobs

aws glue create-database \
    --database-input "{\"Name\":\"${TEST_PREFIX}_db\",\"Description\":\"Test database\"}"

aws glue create-table \
    --database-name ${TEST_PREFIX}_db \
    --table-input '{
        "Name": "test_counts",
        "StorageDescriptor": {
            "Columns": [
                {"Name": "gene_id", "Type": "string"},
                {"Name": "count", "Type": "int"}
            ],
            "Location": "s3://'${TEST_BUCKET}'/raw/counts/",
            "InputFormat": "org.apache.hadoop.mapred.TextInputFormat",
            "OutputFormat": "org.apache.hadoop.hive.ql.io.HiveIgnoreKeyTextOutputFormat",
            "SerdeInfo": {
                "SerializationLibrary": "org.apache.hadoop.hive.serde2.lazy.LazySimpleSerDe"
            }
        }
    }'

aws glue get-table --database-name ${TEST_PREFIX}_db --name test_counts
cd ..
echo "✓ Glue test passed"

# Test 5: Healthcare Standards with S3
echo -e "\n🏥 Test 5: Healthcare Data Standards"
cd healthcare_standards

# Generate FHIR data
python3 << 'EOF'
from fhir_genomics_integration import create_genomic_observation
import json

obs = create_genomic_observation("patient-test", {"BRCA1": 123.4, "TP53": 567.8})
with open("test_fhir.json", "w") as f:
    json.dump(obs, f, indent=2)
print("Generated FHIR observation")
EOF

# Upload to S3
aws s3 cp test_fhir.json s3://${TEST_BUCKET}/fhir/observations/test_fhir.json

# Generate HL7 data
python3 << 'EOF'
from hl7_v2_parser import create_hl7_message
msg = create_hl7_message("patient-test", [{"test": "GLU", "value": 95, "unit": "mg/dL"}])
with open("test_hl7.txt", "w") as f:
    f.write(msg)
print("Generated HL7 message")
EOF

aws s3 cp test_hl7.txt s3://${TEST_BUCKET}/hl7/messages/test_hl7.txt

cd ..
echo "✓ Healthcare standards test passed"

# Summary
echo -e "\n═══════════════════════════════════════"
echo -e "✅ All AWS tests passed!"
echo -e "═══════════════════════════════════════"
echo "Test bucket: s3://${TEST_BUCKET}"
echo "Lambda function: ${TEST_PREFIX}-trigger"
echo "Glue database: ${TEST_PREFIX}_db"
echo -e "═══════════════════════════════════════"
