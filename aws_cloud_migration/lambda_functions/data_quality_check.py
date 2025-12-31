"""
AWS Lambda Function: Data Quality Check

This Lambda function validates RNA-Seq data files for quality issues
and sends SNS notifications on validation failures.
"""

import json
import boto3
import logging
import os
import csv
from io import StringIO
from typing import Dict, Any, List, Tuple, Optional
from botocore.exceptions import ClientError

# Configure logging
logger = logging.getLogger()
logger.setLevel(logging.INFO)

# Initialize AWS clients
s3_client = boto3.client('s3')
sns_client = boto3.client('sns')

# Environment variables
SNS_TOPIC_ARN = os.environ.get('SNS_TOPIC_ARN', '')
MAX_MISSING_PERCENT = float(os.environ.get('MAX_MISSING_PERCENT', '10.0'))
MIN_SAMPLE_COUNT = int(os.environ.get('MIN_SAMPLE_COUNT', '3'))
REQUIRED_COLUMNS = os.environ.get('REQUIRED_COLUMNS', 'gene_id').split(',')


class DataQualityChecker:
    """Class for validating RNA-Seq data quality."""
    
    def __init__(self, bucket: str, key: str):
        """
        Initialize quality checker.
        
        Args:
            bucket: S3 bucket name
            key: S3 object key
        """
        self.bucket = bucket
        self.key = key
        self.errors: List[str] = []
        self.warnings: List[str] = []
        
    def check_file_format(self) -> bool:
        """
        Check if file format is valid (CSV/TSV).
        
        Returns:
            True if valid format, False otherwise
        """
        logger.info(f"Checking file format for: {self.key}")
        
        valid_extensions = ['.csv', '.tsv', '.txt']
        
        if not any(self.key.lower().endswith(ext) for ext in valid_extensions):
            self.errors.append(f"Invalid file format. Expected CSV/TSV, got: {self.key}")
            return False
        
        logger.info("File format check passed")
        return True
    
    def read_file_sample(self, max_lines: int = 100) -> Optional[str]:
        """
        Read first N lines of file from S3.
        
        Args:
            max_lines: Maximum number of lines to read
            
        Returns:
            File content as string, or None if error
        """
        try:
            response = s3_client.get_object(Bucket=self.bucket, Key=self.key)
            
            # Read first N lines
            lines = []
            for i, line in enumerate(response['Body'].iter_lines()):
                if i >= max_lines:
                    break
                lines.append(line.decode('utf-8'))
            
            content = '\n'.join(lines)
            logger.info(f"Read {len(lines)} lines from file")
            
            return content
            
        except ClientError as e:
            error_msg = f"Error reading file from S3: {str(e)}"
            logger.error(error_msg)
            self.errors.append(error_msg)
            return None
    
    def detect_delimiter(self, content: str) -> str:
        """
        Detect delimiter (comma or tab) in file.
        
        Args:
            content: File content
            
        Returns:
            Detected delimiter
        """
        # Try to detect using csv.Sniffer
        try:
            sniffer = csv.Sniffer()
            delimiter = sniffer.sniff(content[:1000]).delimiter
            logger.info(f"Detected delimiter: {repr(delimiter)}")
            return delimiter
        except:
            # Default to comma
            logger.warning("Could not detect delimiter, defaulting to comma")
            return ','
    
    def validate_headers(self, content: str) -> Tuple[bool, List[str]]:
        """
        Validate column headers.
        
        Args:
            content: File content
            
        Returns:
            Tuple of (is_valid, column_names)
        """
        logger.info("Validating column headers")
        
        try:
            delimiter = self.detect_delimiter(content)
            reader = csv.reader(StringIO(content), delimiter=delimiter)
            headers = next(reader)
            
            # Remove whitespace from headers
            headers = [h.strip() for h in headers]
            
            logger.info(f"Found {len(headers)} columns: {headers[:10]}...")
            
            # Check for required columns
            missing_cols = [col for col in REQUIRED_COLUMNS if col not in headers]
            
            if missing_cols:
                error_msg = f"Missing required columns: {missing_cols}"
                logger.error(error_msg)
                self.errors.append(error_msg)
                return False, headers
            
            # Check for duplicate headers
            if len(headers) != len(set(headers)):
                error_msg = "Duplicate column headers found"
                logger.error(error_msg)
                self.errors.append(error_msg)
                return False, headers
            
            # Check for sample columns
            sample_cols = [h for h in headers if h.startswith('sample_') or 
                          h.startswith('Sample') or h not in REQUIRED_COLUMNS]
            
            if len(sample_cols) < MIN_SAMPLE_COUNT:
                error_msg = f"Insufficient sample columns. Found {len(sample_cols)}, minimum required: {MIN_SAMPLE_COUNT}"
                logger.error(error_msg)
                self.errors.append(error_msg)
                return False, headers
            
            logger.info("Header validation passed")
            return True, headers
            
        except Exception as e:
            error_msg = f"Error validating headers: {str(e)}"
            logger.error(error_msg)
            self.errors.append(error_msg)
            return False, []
    
    def check_missing_values(self, content: str, headers: List[str]) -> bool:
        """
        Check for missing values in data.
        
        Args:
            content: File content
            headers: Column headers
            
        Returns:
            True if within acceptable limits, False otherwise
        """
        logger.info("Checking for missing values")
        
        try:
            delimiter = self.detect_delimiter(content)
            reader = csv.DictReader(StringIO(content), delimiter=delimiter)
            
            total_cells = 0
            missing_cells = 0
            row_count = 0
            
            for row in reader:
                row_count += 1
                for col in headers:
                    total_cells += 1
                    value = row.get(col, '').strip()
                    
                    if value == '' or value.lower() in ['na', 'nan', 'null', 'none', 'n/a']:
                        missing_cells += 1
            
            if total_cells == 0:
                error_msg = "No data rows found in file"
                logger.error(error_msg)
                self.errors.append(error_msg)
                return False
            
            missing_percent = (missing_cells / total_cells) * 100
            
            logger.info(f"Data statistics: {row_count} rows, {total_cells} cells, "
                       f"{missing_cells} missing ({missing_percent:.2f}%)")
            
            if missing_percent > MAX_MISSING_PERCENT:
                error_msg = (f"Too many missing values: {missing_percent:.2f}% "
                           f"(threshold: {MAX_MISSING_PERCENT}%)")
                logger.error(error_msg)
                self.errors.append(error_msg)
                return False
            
            if missing_percent > 5.0:
                warning_msg = f"Warning: {missing_percent:.2f}% missing values detected"
                logger.warning(warning_msg)
                self.warnings.append(warning_msg)
            
            logger.info("Missing values check passed")
            return True
            
        except Exception as e:
            error_msg = f"Error checking missing values: {str(e)}"
            logger.error(error_msg)
            self.errors.append(error_msg)
            return False
    
    def validate_numeric_data(self, content: str, headers: List[str]) -> bool:
        """
        Validate that count columns contain numeric data.
        
        Args:
            content: File content
            headers: Column headers
            
        Returns:
            True if valid, False otherwise
        """
        logger.info("Validating numeric data in count columns")
        
        try:
            delimiter = self.detect_delimiter(content)
            reader = csv.DictReader(StringIO(content), delimiter=delimiter)
            
            # Get sample columns (exclude metadata columns)
            sample_cols = [h for h in headers if h not in REQUIRED_COLUMNS]
            
            non_numeric_count = 0
            total_count = 0
            
            for row_num, row in enumerate(reader):
                for col in sample_cols:
                    value = row.get(col, '').strip()
                    total_count += 1
                    
                    if value and value.lower() not in ['na', 'nan', 'null', 'none', 'n/a']:
                        try:
                            float(value)
                        except ValueError:
                            non_numeric_count += 1
                            if non_numeric_count <= 5:  # Log first 5 errors
                                logger.warning(f"Non-numeric value '{value}' in row {row_num + 2}, column {col}")
            
            if non_numeric_count > 0:
                non_numeric_percent = (non_numeric_count / total_count) * 100
                
                if non_numeric_percent > 1.0:  # More than 1% non-numeric
                    error_msg = (f"Too many non-numeric values in count columns: "
                               f"{non_numeric_count} ({non_numeric_percent:.2f}%)")
                    logger.error(error_msg)
                    self.errors.append(error_msg)
                    return False
                else:
                    warning_msg = f"Found {non_numeric_count} non-numeric values in count columns"
                    logger.warning(warning_msg)
                    self.warnings.append(warning_msg)
            
            logger.info("Numeric data validation passed")
            return True
            
        except Exception as e:
            error_msg = f"Error validating numeric data: {str(e)}"
            logger.error(error_msg)
            self.errors.append(error_msg)
            return False
    
    def run_all_checks(self) -> Dict[str, Any]:
        """
        Run all quality checks.
        
        Returns:
            Dictionary with check results
        """
        logger.info("Starting data quality checks")
        
        # Check file format
        if not self.check_file_format():
            return {
                'valid': False,
                'errors': self.errors,
                'warnings': self.warnings
            }
        
        # Read file sample
        content = self.read_file_sample()
        if content is None:
            return {
                'valid': False,
                'errors': self.errors,
                'warnings': self.warnings
            }
        
        # Validate headers
        headers_valid, headers = self.validate_headers(content)
        if not headers_valid:
            return {
                'valid': False,
                'errors': self.errors,
                'warnings': self.warnings
            }
        
        # Check missing values
        if not self.check_missing_values(content, headers):
            return {
                'valid': False,
                'errors': self.errors,
                'warnings': self.warnings
            }
        
        # Validate numeric data
        if not self.validate_numeric_data(content, headers):
            return {
                'valid': False,
                'errors': self.errors,
                'warnings': self.warnings
            }
        
        logger.info("All quality checks passed")
        return {
            'valid': True,
            'errors': self.errors,
            'warnings': self.warnings
        }


def send_sns_notification(
    topic_arn: str,
    subject: str,
    message: str
) -> bool:
    """
    Send SNS notification.
    
    Args:
        topic_arn: SNS topic ARN
        subject: Email subject
        message: Email message
        
    Returns:
        True if successful, False otherwise
    """
    if not topic_arn:
        logger.warning("SNS topic ARN not configured, skipping notification")
        return False
    
    try:
        logger.info(f"Sending SNS notification to: {topic_arn}")
        
        response = sns_client.publish(
            TopicArn=topic_arn,
            Subject=subject,
            Message=message
        )
        
        logger.info(f"SNS notification sent - Message ID: {response['MessageId']}")
        return True
        
    except ClientError as e:
        logger.error(f"Error sending SNS notification: {str(e)}")
        return False


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
        # Extract S3 information
        if 'Records' in event and len(event['Records']) > 0:
            record = event['Records'][0]
            bucket = record['s3']['bucket']['name']
            key = record['s3']['object']['key']
        else:
            # Direct invocation format
            bucket = event.get('bucket', '')
            key = event.get('key', '')
        
        if not bucket or not key:
            return {
                'statusCode': 400,
                'body': json.dumps({
                    'message': 'Invalid event structure',
                    'error': 'Missing bucket or key'
                })
            }
        
        logger.info(f"Validating file: s3://{bucket}/{key}")
        
        # Run quality checks
        checker = DataQualityChecker(bucket, key)
        results = checker.run_all_checks()
        
        # Send notification if validation failed
        if not results['valid']:
            subject = f"Data Quality Check Failed: {key}"
            message = f"""
Data Quality Check Failed

File: s3://{bucket}/{key}

Errors:
{chr(10).join('- ' + e for e in results['errors'])}

Warnings:
{chr(10).join('- ' + w for w in results['warnings'])}

Please review the file and re-upload if necessary.
            """.strip()
            
            send_sns_notification(SNS_TOPIC_ARN, subject, message)
            
            return {
                'statusCode': 422,
                'body': json.dumps({
                    'message': 'Data quality check failed',
                    'bucket': bucket,
                    'key': key,
                    'results': results
                })
            }
        
        # Send notification for warnings
        if results['warnings']:
            subject = f"Data Quality Warnings: {key}"
            message = f"""
Data Quality Check Passed with Warnings

File: s3://{bucket}/{key}

Warnings:
{chr(10).join('- ' + w for w in results['warnings'])}

The file passed validation but you may want to review these warnings.
            """.strip()
            
            send_sns_notification(SNS_TOPIC_ARN, subject, message)
        
        # Success response
        return {
            'statusCode': 200,
            'body': json.dumps({
                'message': 'Data quality check passed',
                'bucket': bucket,
                'key': key,
                'results': results
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
    test_event = {
        "bucket": "bioinformatics-pipeline",
        "key": "raw/counts/sample_counts.csv"
    }
    
    class MockContext:
        def __init__(self):
            self.function_name = "data_quality_check"
    
    result = lambda_handler(test_event, MockContext())
    print(json.dumps(result, indent=2))
