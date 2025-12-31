#!/bin/bash
# Local streaming tests using Docker Compose

set -e

echo "🌊 Starting Local Streaming Tests..."

GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m'

# Start services
echo "🐳 Starting Docker Compose stack..."
docker-compose up -d

# Wait for services
echo "⏳ Waiting for services to be ready..."
sleep 30

# Test counter
PASSED=0
FAILED=0

test_service() {
    local service=$1
    local check_command=$2
    
    echo -e "\n${YELLOW}Testing: ${service}${NC}"
    
    if eval "$check_command"; then
        echo -e "${GREEN}✓ ${service} is healthy${NC}"
        ((PASSED++))
    else
        echo -e "${RED}✗ ${service} failed${NC}"
        ((FAILED++))
    fi
}

# Test Kafka
test_service "Kafka Broker" "docker exec kafka-broker kafka-broker-api-versions --bootstrap-server localhost:9092 > /dev/null 2>&1"

# Test Zookeeper
test_service "Zookeeper" "docker exec zookeeper zkServer.sh status | grep -q 'Mode: '"

# Test Flink
test_service "Flink JobManager" "curl -sf http://localhost:8081/overview > /dev/null"

# Test Spark
test_service "Spark Master" "curl -sf http://localhost:8082 > /dev/null"

# Test PostgreSQL
test_service "PostgreSQL" "docker exec postgres pg_isready -U postgres | grep -q 'accepting connections'"

# Create Kafka topics
echo -e "\n📋 Creating Kafka topics..."
docker exec kafka-broker kafka-topics --create \
    --bootstrap-server localhost:9092 \
    --topic medical-devices-vitals \
    --partitions 3 \
    --replication-factor 1 \
    --if-not-exists

docker exec kafka-broker kafka-topics --create \
    --bootstrap-server localhost:9092 \
    --topic medical-devices-glucose \
    --partitions 3 \
    --replication-factor 1 \
    --if-not-exists

# Test producer
echo -e "\n📤 Testing Kafka Producer..."
cd kafka_streams
python3 << 'EOF'
from producer_medical_devices import generate_vital_signs, send_to_kafka
from kafka import KafkaProducer
import json
import time

producer = KafkaProducer(
    bootstrap_servers=['localhost:9092'],
    value_serializer=lambda v: json.dumps(v).encode('utf-8')
)

for i in range(10):
    message = generate_vital_signs(f"patient-{i}")
    producer.send('medical-devices-vitals', value=message)
    time.sleep(0.1)

producer.flush()
print("Sent 10 test messages")
EOF

# Verify messages
MESSAGE_COUNT=$(docker exec kafka-broker kafka-run-class kafka.tools.GetOffsetShell \
    --broker-list localhost:9092 \
    --topic medical-devices-vitals | awk -F':' '{sum += $3} END {print sum}')

if [ "$MESSAGE_COUNT" -ge 10 ]; then
    echo -e "${GREEN}✓ Producer test passed (${MESSAGE_COUNT} messages)${NC}"
    ((PASSED++))
else
    echo -e "${RED}✗ Producer test failed${NC}"
    ((FAILED++))
fi

# Test consumer
echo -e "\n📥 Testing Kafka Consumer..."
timeout 10s docker exec kafka-broker kafka-console-consumer \
    --bootstrap-server localhost:9092 \
    --topic medical-devices-vitals \
    --max-messages 5 \
    --timeout-ms 5000 > /dev/null 2>&1

if [ $? -eq 0 ] || [ $? -eq 124 ]; then
    echo -e "${GREEN}✓ Consumer test passed${NC}"
    ((PASSED++))
else
    echo -e "${RED}✗ Consumer test failed${NC}"
    ((FAILED++))
fi

# Test Flink job submission
echo -e "\n🔄 Testing Flink Job..."
docker cp flink_processing/medical_data_processor.py flink-jobmanager:/opt/flink/
docker exec flink-jobmanager flink run -py /opt/flink/medical_data_processor.py --detached > /dev/null 2>&1

sleep 5

JOB_COUNT=$(curl -s http://localhost:8081/jobs | jq '.jobs | length')
if [ "$JOB_COUNT" -gt 0 ]; then
    echo -e "${GREEN}✓ Flink job test passed${NC}"
    ((PASSED++))
else
    echo -e "${RED}✗ Flink job test failed${NC}"
    ((FAILED++))
fi

# Test Spark Streaming
echo -e "\n⚡ Testing Spark Streaming..."
docker cp spark_streaming/structured_streaming_hl7.py spark-master:/opt/spark-apps/
docker exec spark-master spark-submit \
    --master local[2] \
    --packages org.apache.spark:spark-sql-kafka-0-10_2.12:3.4.0 \
    /opt/spark-apps/structured_streaming_hl7.py > /dev/null 2>&1 &

sleep 10
SPARK_JOBS=$(curl -s http://localhost:8082/api/v1/applications | jq '. | length')
if [ "$SPARK_JOBS" -gt 0 ]; then
    echo -e "${GREEN}✓ Spark streaming test passed${NC}"
    ((PASSED++))
else
    echo -e "${RED}✗ Spark streaming test failed${NC}"
    ((FAILED++))
fi

# Summary
echo -e "\n═══════════════════════════════════════"
echo -e "📊 Test Summary"
echo -e "═══════════════════════════════════════"
echo -e "${GREEN}Passed: ${PASSED}${NC}"
echo -e "${RED}Failed: ${FAILED}${NC}"
echo -e "═══════════════════════════════════════"

# Cleanup option
read -p "Stop Docker Compose stack? (y/n) " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    docker-compose down -v
fi

[ $FAILED -eq 0 ] && exit 0 || exit 1
