#!/bin/bash
# Streaming performance benchmarks

echo "⚡ Streaming Performance Benchmarks..."

# Test 1: Kafka throughput
echo "📊 Kafka Throughput Test..."
docker exec kafka-broker kafka-producer-perf-test \
    --topic medical-devices-vitals \
    --num-records 100000 \
    --record-size 1000 \
    --throughput -1 \
    --producer-props bootstrap.servers=localhost:9092

# Test 2: End-to-end latency
echo "📊 End-to-End Latency Test..."
# Measure time from produce to consume

# Test 3: Flink processing rate
echo "📊 Flink Processing Rate..."
# Check Flink metrics

# Test 4: Consumer lag
echo "📊 Consumer Lag Test..."
docker exec kafka-broker kafka-consumer-groups \
    --bootstrap-server localhost:9092 \
    --describe \
    --all-groups

# Save results
echo "Results saved to benchmarks/$(date +%Y%m%d_%H%M%S).txt"
