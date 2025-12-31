# Comprehensive Testing Infrastructure

This repository now includes complete automated testing infrastructure for AWS cloud migration and healthcare streaming components.

## 📁 Structure

### 1. AWS Cloud Migration (`aws_cloud_migration/`)

Comprehensive testing infrastructure for AWS-based bioinformatics pipelines:

- **Healthcare Standards**: FHIR, HL7 v2, and OMOP CDM implementations
- **AWS Services**: Lambda, Glue, S3, ECS, EMR integration
- **Test Scripts**: Local (no AWS) and AWS integration tests
- **CI/CD**: GitHub Actions workflow for automated testing

#### Quick Start - AWS Cloud Migration
```bash
# Local testing (no AWS account required)
cd aws_cloud_migration
bash scripts/run_local_tests.sh

# AWS integration testing (requires credentials)
aws configure
bash scripts/run_aws_tests.sh

# Performance benchmarks
bash scripts/performance_benchmark.sh
```

#### What's Tested
- ✅ Healthcare data standards (FHIR, HL7, OMOP)
- ✅ Docker container builds
- ✅ Terraform configuration validation
- ✅ Python code quality and syntax
- ✅ AWS Lambda functions (mocked)
- ✅ AWS Glue jobs
- ✅ S3 data lake operations

### 2. Healthcare Streaming (`enterprise-ai-workflows/project4-healthcare-streaming/`)

Real-time streaming infrastructure for healthcare data:

- **Streaming Platforms**: Kafka, Apache Flink, Spark Structured Streaming
- **Data Generators**: Medical device simulators
- **Test Scripts**: Docker Compose-based local testing
- **CI/CD**: GitHub Actions workflow with integration tests

#### Quick Start - Healthcare Streaming
```bash
# Local streaming tests
cd enterprise-ai-workflows/project4-healthcare-streaming
bash scripts/run_local_tests.sh

# Performance benchmarks
bash scripts/performance_benchmark.sh
```

#### What's Tested
- ✅ Kafka broker health
- ✅ Zookeeper coordination
- ✅ Flink job submission
- ✅ Spark streaming jobs
- ✅ PostgreSQL database
- ✅ Producer/consumer functionality
- ✅ End-to-end data flow

## 🔄 GitHub Actions CI/CD

### AWS Cloud Migration Tests (`.github/workflows/aws-cloud-tests.yml`)
Triggered on push/PR to `aws_cloud_migration/**`:
- Local tests (no AWS credentials)
- Terraform validation
- Docker builds
- AWS integration tests (main branch only)

### Streaming Tests (`.github/workflows/streaming-tests.yml`)
Triggered on push/PR to `project4-healthcare-streaming/**`:
- Code quality checks (flake8, black, isort)
- Unit tests with coverage
- Docker Compose integration tests
- AWS deployment validation (main branch only)

## 📊 Test Coverage

### AWS Cloud Migration
- **Unit Tests**: Healthcare standards modules
- **Integration Tests**: AWS service interactions
- **Validation Tests**: Terraform, Docker, JSON configs
- **Mock Tests**: Lambda functions without AWS

### Healthcare Streaming
- **Service Health**: All containers start successfully
- **Message Flow**: Producer → Kafka → Consumer
- **Processing**: Flink and Spark job execution
- **Performance**: Throughput and latency benchmarks

## 🚀 Running Tests Locally

### Prerequisites
```bash
# Python dependencies
pip install -r aws_cloud_migration/tests/requirements.txt
pip install -r enterprise-ai-workflows/project4-healthcare-streaming/tests/requirements.txt

# For Docker tests
docker --version

# For Terraform tests
terraform --version

# For AWS integration tests
aws configure
```

### Test Commands

**AWS Cloud Migration - Full Suite:**
```bash
cd aws_cloud_migration
python3 -m pytest tests/ -v --cov=. --cov-report=term-missing
```

**Healthcare Streaming - With Docker:**
```bash
cd enterprise-ai-workflows/project4-healthcare-streaming
docker-compose up -d
bash scripts/run_local_tests.sh
docker-compose down -v
```

## 📝 Test Documentation

- `aws_cloud_migration/TESTING.md` - AWS cloud testing guide
- `enterprise-ai-workflows/project4-healthcare-streaming/TESTING.md` - Streaming testing guide

## 🎯 Success Criteria

✅ All test scripts are executable  
✅ Tests run successfully in CI/CD  
✅ Clear documentation for running tests  
✅ Both local and AWS testing covered  
✅ Performance benchmarking included  
✅ Proper error handling and cleanup  
✅ Test results logged and reportable  

## 🔧 Maintenance

### Adding New Tests

**AWS Cloud Migration:**
1. Add test functions to `aws_cloud_migration/tests/test_*.py`
2. Update `run_local_tests.sh` for new test scenarios
3. Ensure CI workflow catches new test patterns

**Healthcare Streaming:**
1. Add tests to `project4-healthcare-streaming/tests/`
2. Update `run_local_tests.sh` for integration scenarios
3. Add new services to `docker-compose.yml` if needed

## 🐛 Troubleshooting

**Common Issues:**

1. **Docker permission denied**: Add user to docker group
   ```bash
   sudo usermod -aG docker $USER
   ```

2. **AWS credentials not found**: Configure AWS CLI
   ```bash
   aws configure
   ```

3. **Terraform backend errors**: Use `-backend=false` for validation
   ```bash
   terraform init -backend=false
   ```

4. **Port conflicts in Docker**: Stop conflicting containers
   ```bash
   docker ps
   docker stop <container-id>
   ```

## 📈 Future Enhancements

- [ ] Add more comprehensive integration tests
- [ ] Implement end-to-end workflow tests
- [ ] Add performance regression testing
- [ ] Create test data generators
- [ ] Implement chaos engineering tests
- [ ] Add security scanning to CI/CD

## 📄 License

MIT License - See main repository LICENSE file
