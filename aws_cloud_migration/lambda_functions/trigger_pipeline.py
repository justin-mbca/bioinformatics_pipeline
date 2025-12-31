"""
AWS Lambda Function - Pipeline Trigger
Triggered by S3 events to start Glue jobs
"""
import json
import boto3
import os


def lambda_handler(event, context):
    """
    Lambda handler to trigger Glue job when new data arrives in S3
    
    Args:
        event: S3 event notification
        context: Lambda context
        
    Returns:
        Status code and message
    """
    glue = boto3.client('glue')
    
    try:
        # Parse S3 event
        for record in event['Records']:
            bucket = record['s3']['bucket']['name']
            key = record['s3']['object']['key']
            
            print(f"Processing S3 object: s3://{bucket}/{key}")
            
            # Start Glue job
            response = glue.start_job_run(
                JobName=os.environ.get('GLUE_JOB_NAME', 'bioinformatics-etl'),
                Arguments={
                    '--S3_INPUT_PATH': f"s3://{bucket}/{key}",
                    '--S3_OUTPUT_PATH': f"s3://{bucket}/processed/"
                }
            )
            
            job_run_id = response['JobRunId']
            print(f"Started Glue job run: {job_run_id}")
            
        return {
            'statusCode': 200,
            'body': json.dumps(f'Successfully triggered pipeline')
        }
        
    except Exception as e:
        print(f"Error: {str(e)}")
        return {
            'statusCode': 500,
            'body': json.dumps(f'Error: {str(e)}')
        }


if __name__ == "__main__":
    # Test locally
    test_event = {
        "Records": [{
            "s3": {
                "bucket": {"name": "test-bucket"},
                "object": {"key": "raw/test.csv"}
            }
        }]
    }
    result = lambda_handler(test_event, None)
    print(result)
