#!/bin/bash
#
# EMR Bootstrap Actions for Bioinformatics Pipeline
# 
# This script installs Python packages and configures Spark settings
# for RNA-Seq analysis on EMR clusters.
#

set -e  # Exit on error
set -x  # Print commands

echo "========================================"
echo "Starting EMR Bootstrap Actions"
echo "========================================"

# Update system packages
echo "Updating system packages..."
sudo yum update -y

# Install Python development tools
echo "Installing Python development tools..."
sudo yum install -y python3-devel gcc gcc-c++ make

# Upgrade pip
echo "Upgrading pip..."
sudo python3 -m pip install --upgrade pip setuptools wheel

# Install bioinformatics Python packages
echo "Installing bioinformatics packages..."
sudo python3 -m pip install \
    numpy==1.24.3 \
    pandas==2.0.2 \
    scipy==1.10.1 \
    scikit-learn==1.2.2 \
    biopython==1.81 \
    pysam==0.21.0 \
    htseq==2.0.2 \
    statsmodels==0.14.0

# Install AWS SDK for Python
echo "Installing AWS SDK..."
sudo python3 -m pip install \
    boto3==1.26.137 \
    botocore==1.29.137

# Install additional data science packages
echo "Installing data science packages..."
sudo python3 -m pip install \
    matplotlib==3.7.1 \
    seaborn==0.12.2 \
    plotly==5.14.1

# Configure Spark settings
echo "Configuring Spark settings..."

# Create Spark configuration directory if it doesn't exist
sudo mkdir -p /etc/spark/conf.dist

# Set Spark defaults
cat << 'EOF' | sudo tee -a /etc/spark/conf.dist/spark-defaults.conf

# Bioinformatics Pipeline Spark Configuration
spark.python.profile true
spark.python.profile.memory true
spark.python.worker.memory 4g
spark.executor.memoryOverhead 2048
spark.driver.memoryOverhead 2048

# Optimize for genomics data processing
spark.sql.files.maxPartitionBytes 256m
spark.sql.shuffle.partitions 200
spark.default.parallelism 200

# Enable advanced optimizations
spark.sql.adaptive.enabled true
spark.sql.adaptive.coalescePartitions.enabled true
spark.sql.adaptive.skewJoin.enabled true

# Network and timeout settings for large data transfers
spark.network.timeout 800s
spark.executor.heartbeatInterval 60s
spark.storage.blockManagerSlaveTimeoutMs 800000

# Serialization
spark.serializer org.apache.spark.serializer.KryoSerializer
spark.kryoserializer.buffer.max 512m

# History server
spark.eventLog.enabled true
spark.eventLog.dir hdfs:///var/log/spark/apps
spark.history.fs.logDirectory hdfs:///var/log/spark/apps
EOF

# Configure Python environment
echo "Configuring Python environment..."
cat << 'EOF' | sudo tee -a /etc/profile.d/bioinformatics.sh

# Bioinformatics Pipeline Environment Variables
export PYTHONPATH=/usr/lib/spark/python:$PYTHONPATH
export PYSPARK_PYTHON=/usr/bin/python3
export PYSPARK_DRIVER_PYTHON=/usr/bin/python3

# Increase file descriptor limits for genomics data
ulimit -n 65536
EOF

# Source the environment
source /etc/profile.d/bioinformatics.sh

# Create working directories
echo "Creating working directories..."
sudo mkdir -p /mnt/tmp
sudo chmod 777 /mnt/tmp
sudo mkdir -p /mnt/data
sudo chmod 777 /mnt/data

# Set up logging
echo "Configuring logging..."
sudo mkdir -p /var/log/bioinformatics
sudo chmod 777 /var/log/bioinformatics

# Configure log4j for better logging
cat << 'EOF' | sudo tee /etc/spark/conf.dist/log4j.properties
# Root logger option
log4j.rootCategory=INFO, console, file

# Console appender
log4j.appender.console=org.apache.log4j.ConsoleAppender
log4j.appender.console.target=System.err
log4j.appender.console.layout=org.apache.log4j.PatternLayout
log4j.appender.console.layout.ConversionPattern=%d{yy/MM/dd HH:mm:ss} %p %c{1}: %m%n

# File appender
log4j.appender.file=org.apache.log4j.RollingFileAppender
log4j.appender.file.File=/var/log/bioinformatics/spark.log
log4j.appender.file.MaxFileSize=100MB
log4j.appender.file.MaxBackupIndex=10
log4j.appender.file.layout=org.apache.log4j.PatternLayout
log4j.appender.file.layout.ConversionPattern=%d{yy/MM/dd HH:mm:ss} %p %c{1}: %m%n

# Reduce verbosity for specific loggers
log4j.logger.org.apache.spark.repl.Main=WARN
log4j.logger.org.spark_project.jetty=WARN
log4j.logger.org.spark_project.jetty.util.component.AbstractLifeCycle=ERROR
log4j.logger.org.apache.spark.repl.SparkIMain$exprTyper=INFO
log4j.logger.org.apache.spark.repl.SparkILoop$SparkILoopInterpreter=INFO
log4j.logger.org.apache.parquet=ERROR
log4j.logger.parquet=ERROR
log4j.logger.org.apache.hadoop.hive.metastore.RetryingHMSHandler=FATAL
log4j.logger.org.apache.hadoop.hive.ql.exec.FunctionRegistry=ERROR
EOF

# Install system monitoring tools
echo "Installing monitoring tools..."
sudo yum install -y htop iotop sysstat

# Download analysis scripts from S3 (if available)
echo "Downloading analysis scripts..."
if aws s3 ls s3://bioinformatics-pipeline/emr-scripts/ 2>/dev/null; then
    sudo mkdir -p /opt/bioinformatics
    aws s3 sync s3://bioinformatics-pipeline/emr-scripts/ /opt/bioinformatics/
    sudo chmod +x /opt/bioinformatics/*.py
    sudo chmod +x /opt/bioinformatics/*.sh
    echo "Analysis scripts downloaded successfully"
else
    echo "No scripts found in S3, skipping download"
fi

# Verify installations
echo "Verifying installations..."
python3 --version
pip3 --version
python3 -c "import numpy; print(f'NumPy version: {numpy.__version__}')"
python3 -c "import pandas; print(f'Pandas version: {pandas.__version__}')"
python3 -c "import scipy; print(f'SciPy version: {scipy.__version__}')"
python3 -c "import Bio; print(f'Biopython version: {Bio.__version__}')"
python3 -c "import boto3; print(f'Boto3 version: {boto3.__version__}')"

echo "========================================"
echo "Bootstrap Actions Completed Successfully"
echo "========================================"
