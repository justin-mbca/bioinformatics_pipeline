#!/bin/bash
# Local testing script - no AWS account needed
# Tests: Python scripts, Docker containers, Terraform validation, Healthcare standards

set -e  # Exit on error

echo "🧪 Starting Local Tests for AWS Cloud Migration Tools..."

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Test counter
TESTS_PASSED=0
TESTS_FAILED=0

# Function to run test
run_test() {
    local test_name="$1"
    local test_command="$2"
    
    echo -e "\n${YELLOW}Testing: ${test_name}${NC}"
    
    if eval "$test_command"; then
        echo -e "${GREEN}✓ PASSED${NC}"
        ((TESTS_PASSED++))
        return 0
    else
        echo -e "${RED}✗ FAILED${NC}"
        ((TESTS_FAILED++))
        return 1
    fi
}

# Setup virtual environment
echo "Setting up Python environment..."
python3 -m venv .venv
source .venv/bin/activate
pip install -q boto3 pandas numpy biopython pytest pytest-mock pytest-cov moto

# Test 1: Python syntax validation
run_test "Python Syntax Check" "python3 -m py_compile glue_jobs/*.py lambda_functions/*.py healthcare_standards/*.py s3_data_lake/*.py"

# Test 2: Healthcare Standards - FHIR
run_test "FHIR Genomics Integration" "python3 -c 'from healthcare_standards.fhir_genomics_integration import create_genomic_observation; result = create_genomic_observation(\"patient-123\", {\"BRCA1\": 234.5}); assert result is not None'"

# Test 3: Healthcare Standards - HL7
run_test "HL7 v2 Parser" "python3 -c 'from healthcare_standards.hl7_v2_parser import parse_hl7_message; msg = \"MSH|^~\\&|LAB|HOSP|||20231231120000||ORU^R01|001|P|2.5\"; result = parse_hl7_message(msg); assert result is not None'"

# Test 4: Healthcare Standards - OMOP
run_test "OMOP CDM Mapping" "python3 -c 'from healthcare_standards.omop_cdm_mapping import map_to_omop; result = map_to_omop({\"patient_id\": \"123\"}); assert result is not None'"

# Test 5: Unit tests
run_test "Unit Tests" "pytest tests/ -v --cov=. --cov-report=term-missing"

# Test 6: Docker build
run_test "Docker Container Build" "docker build -t rnaseq-test:latest -f ecs_containers/Dockerfile.rnaseq ."

# Test 7: Docker container run
run_test "Docker Container Run" "docker run --rm -e TEST_MODE=true rnaseq-test:latest python -c 'print(\"Container works\")'"

# Test 8: Terraform validation
run_test "Terraform Syntax" "cd terraform && terraform init -backend=false && terraform validate"

# Test 9: Terraform format check
run_test "Terraform Formatting" "cd terraform && terraform fmt -check -recursive"

# Test 10: JSON validation
run_test "JSON Configuration Files" "python3 -m json.tool glue_jobs/glue_job_definition.json > /dev/null && python3 -m json.tool emr_spark/emr_cluster_config.json > /dev/null"

# Test 11: Lambda function mock test
run_test "Lambda Trigger Function (Mock)" "python3 -c '
import sys
sys.path.insert(0, \"lambda_functions\")
from unittest.mock import MagicMock, patch
import trigger_pipeline

with patch(\"trigger_pipeline.boto3.client\") as mock_client:
    mock_glue = MagicMock()
    mock_glue.start_job_run.return_value = {\"JobRunId\": \"test-123\"}
    mock_client.return_value = mock_glue
    
    event = {\"Records\": [{\"s3\": {\"bucket\": {\"name\": \"test\"}, \"object\": {\"key\": \"test.csv\"}}}]}
    result = trigger_pipeline.lambda_handler(event, None)
    assert result[\"statusCode\"] == 200
'"

# Test 12: Data quality checks
run_test "Sample Data Validation" "python3 -c '
import pandas as pd
import os
if os.path.exists(\"../../data/sample_counts.csv\"):
    df = pd.read_csv(\"../../data/sample_counts.csv\")
    assert len(df) > 0, \"Sample data is empty\"
    assert not df.isnull().all().any(), \"All null columns found\"
    print(f\"Sample data has {len(df)} rows, {len(df.columns)} columns\")
else:
    print(\"No sample data found - skipping\")
'"

# Cleanup
deactivate
docker rmi rnaseq-test:latest 2>/dev/null || true

# Summary
echo -e "\n═══════════════════════════════════════"
echo -e "📊 Test Summary"
echo -e "═══════════════════════════════════════"
echo -e "${GREEN}Passed: ${TESTS_PASSED}${NC}"
echo -e "${RED}Failed: ${TESTS_FAILED}${NC}"
echo -e "═══════════════════════════════════════"

if [ $TESTS_FAILED -eq 0 ]; then
    echo -e "${GREEN}✓ All tests passed!${NC}"
    exit 0
else
    echo -e "${RED}✗ Some tests failed${NC}"
    exit 1
fi
