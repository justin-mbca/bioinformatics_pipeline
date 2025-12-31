"""
AWS Lambda Function: S3 Event-Triggered Pipeline

This Lambda function is triggered by S3 events when new RNA-Seq files are uploaded.
It starts the AWS Glue ETL job to process the new data.
"""

import json
import boto3
import logging
import os
from typing import Dict, Any, Optional
from botocore.exceptions import ClientError
from datetime import datetime

# Configure logging
logger = logging.getLogger()
logger.setLevel(logging.INFO)

# Initialize AWS clients
glue_client = boto3.client('glue')
s3_client = boto3.client('s3')

# Environment variables
GLUE_JOB_NAME = os.environ.get('GLUE_JOB_NAME', 'rnaseq-etl-pipeline')
GLUE_DATABASE = os.environ.get('GLUE_DATABASE', 'rnaseq_db')
GLUE_TABLE = os.environ.get('GLUE_TABLE', 'raw_counts')
OUTPUT_BUCKET = os.environ.get('OUTPUT_BUCKET', 'bioinformatics-pipeline')
COUNT_THRESHOLD = os.environ.get('COUNT_THRESHOLD', '10')


def extract_s3_info(event: Dict[str, Any]) -> Optional[Dict[str, str]]:
    """
    Extract S3 bucket and key from Lambda event.
    
    Args:
        event: Lambda event dictionary
        
    Returns:
        Dictionary with bucket and key, or None if invalid
    """
    try:
        # S3 event structure
        if 'Records' in event and len(event['Records']) > 0:
            record = event['Records'][0]
            
            if 's3' in record:
                bucket = record['s3']['bucket']['name']
                key = record['s3']['object']['key']
                
                logger.info(f"Extracted S3 info - Bucket: {bucket}, Key: {key}")
                return {'bucket': bucket, 'key': key}
        
        logger.error("Invalid event structure - no S3 records found")
        return None
        
    except (KeyError, IndexError) as e:
        logger.error(f"Error extracting S3 info from event: {str(e)}")
        return None


def validate_file(bucket: str, key: str) -> bool:
    """
    Validate that the uploaded file is a valid RNA-Seq data file.
    
    Args:
        bucket: S3 bucket name
        key: S3 object key
        
    Returns:
        True if valid, False otherwise
    """
    try:
        # Check file extension
        valid_extensions = ['.csv', '.tsv', '.txt', '.parquet']
        if not any(key.lower().endswith(ext) for ext in valid_extensions):
            logger.warning(f"Invalid file extension for key: {key}")
            return False
        
        # Check file size (must be > 0)
        response = s3_client.head_object(Bucket=bucket, Key=key)
        file_size = response['ContentLength']
        
        if file_size == 0:
            logger.warning(f"Empty file detected: {key}")
            return False
        
        logger.info(f"File validation passed - Size: {file_size} bytes")
        return True
        
    except ClientError as e:
        logger.error(f"Error validating file: {str(e)}")
        return False


def start_glue_job(
    bucket: str,
    key: str,
    job_name: str,
    database: str,
    table: str,
    count_threshold: str
) -> Optional[str]:
    """
    Start AWS Glue ETL job with parameters.
    
    Args:
        bucket: S3 bucket containing input data
        key: S3 object key
        job_name: Name of the Glue job to start
        database: Glue database name
        table: Glue table name
        count_threshold: Count threshold for filtering
        
    Returns:
        Job run ID if successful, None otherwise
    """
    try:
        # Construct output path with timestamp
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        output_path = f"s3://{OUTPUT_BUCKET}/processed/normalized/{timestamp}/"
        
        # Job arguments
        arguments = {
            '--database_name': database,
            '--table_name': table,
            '--output_path': output_path,
            '--count_threshold': count_threshold,
            '--input_bucket': bucket,
            '--input_key': key,
            '--job-bookmark-option': 'job-bookmark-enable',
            '--enable-metrics': 'true',
            '--enable-continuous-cloudwatch-log': 'true'
        }
        
        logger.info(f"Starting Glue job: {job_name}")
        logger.info(f"Job arguments: {json.dumps(arguments, indent=2)}")
        
        # Start the job
        response = glue_client.start_job_run(
            JobName=job_name,
            Arguments=arguments,
            Timeout=2880,  # 48 hours
            MaxCapacity=10.0,
            NotificationProperty={
                'NotifyDelayAfter': 60  # Notify after 60 minutes
            }
        )
        
        job_run_id = response['JobRunId']
        logger.info(f"Successfully started Glue job - Run ID: {job_run_id}")
        
        return job_run_id
        
    except ClientError as e:
        error_code = e.response['Error']['Code']
        error_message = e.response['Error']['Message']
        
        if error_code == 'ConcurrentRunsExceededException':
            logger.warning(f"Max concurrent runs exceeded for job: {job_name}")
        else:
            logger.error(f"Error starting Glue job: {error_code} - {error_message}")
        
        return None
    
    except Exception as e:
        logger.error(f"Unexpected error starting Glue job: {str(e)}")
        return None


def get_job_status(job_name: str, job_run_id: str) -> Optional[str]:
    """
    Get current status of a Glue job run.
    
    Args:
        job_name: Name of the Glue job
        job_run_id: Job run ID
        
    Returns:
        Job status string or None if error
    """
    try:
        response = glue_client.get_job_run(
            JobName=job_name,
            RunId=job_run_id
        )
        
        status = response['JobRun']['JobRunState']
        logger.info(f"Job {job_run_id} status: {status}")
        
        return status
        
    except ClientError as e:
        logger.error(f"Error getting job status: {str(e)}")
        return None


def lambda_handler(event: Dict[str, Any], context: Any) -> Dict[str, Any]:
    """
    Main Lambda handler function.
    
    Args:
        event: Lambda event dictionary
        context: Lambda context object
        
    Returns:
        Response dictionary with status and details
    """
    logger.info(f"Received event: {json.dumps(event)}")
    
    try:
        # Extract S3 information from event
        s3_info = extract_s3_info(event)
        
        if not s3_info:
            return {
                'statusCode': 400,
                'body': json.dumps({
                    'message': 'Invalid event structure',
                    'error': 'Could not extract S3 bucket and key from event'
                })
            }
        
        bucket = s3_info['bucket']
        key = s3_info['key']
        
        # Validate file
        if not validate_file(bucket, key):
            return {
                'statusCode': 400,
                'body': json.dumps({
                    'message': 'File validation failed',
                    'bucket': bucket,
                    'key': key
                })
            }
        
        # Start Glue job
        job_run_id = start_glue_job(
            bucket=bucket,
            key=key,
            job_name=GLUE_JOB_NAME,
            database=GLUE_DATABASE,
            table=GLUE_TABLE,
            count_threshold=COUNT_THRESHOLD
        )
        
        if not job_run_id:
            return {
                'statusCode': 500,
                'body': json.dumps({
                    'message': 'Failed to start Glue job',
                    'bucket': bucket,
                    'key': key
                })
            }
        
        # Get initial job status
        status = get_job_status(GLUE_JOB_NAME, job_run_id)
        
        # Success response
        return {
            'statusCode': 200,
            'body': json.dumps({
                'message': 'Successfully triggered ETL pipeline',
                'bucket': bucket,
                'key': key,
                'glue_job_name': GLUE_JOB_NAME,
                'job_run_id': job_run_id,
                'job_status': status,
                'timestamp': datetime.now().isoformat()
            })
        }
        
    except Exception as e:
        logger.error(f"Unexpected error in Lambda handler: {str(e)}")
        
        return {
            'statusCode': 500,
            'body': json.dumps({
                'message': 'Internal server error',
                'error': str(e)
            })
        }


# For local testing
if __name__ == "__main__":
    # Sample S3 event for testing
    test_event = {
        "Records": [
            {
                "eventVersion": "2.1",
                "eventSource": "aws:s3",
                "eventName": "ObjectCreated:Put",
                "s3": {
                    "bucket": {
                        "name": "bioinformatics-pipeline"
                    },
                    "object": {
                        "key": "raw/counts/sample_counts.csv"
                    }
                }
            }
        ]
    }
    
    # Mock context
    class MockContext:
        def __init__(self):
            self.function_name = "trigger_pipeline"
            self.memory_limit_in_mb = 128
            self.invoked_function_arn = "arn:aws:lambda:us-east-1:123456789012:function:trigger_pipeline"
    
    result = lambda_handler(test_event, MockContext())
    print(json.dumps(result, indent=2))
