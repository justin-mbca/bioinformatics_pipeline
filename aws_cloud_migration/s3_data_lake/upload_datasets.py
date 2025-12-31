"""
S3 Data Lake Upload Script
Upload datasets to S3 with proper organization
"""
import boto3
import argparse
import os
from pathlib import Path


def upload_to_s3(bucket, file_path, s3_key):
    """
    Upload file to S3
    
    Args:
        bucket: S3 bucket name
        file_path: Local file path
        s3_key: S3 object key
    """
    s3_client = boto3.client('s3')
    
    try:
        print(f"Uploading {file_path} to s3://{bucket}/{s3_key}")
        s3_client.upload_file(file_path, bucket, s3_key)
        print(f"✓ Upload successful")
        return True
    except Exception as e:
        print(f"✗ Upload failed: {str(e)}")
        return False


def main():
    parser = argparse.ArgumentParser(description='Upload datasets to S3')
    parser.add_argument('--bucket', required=True, help='S3 bucket name')
    parser.add_argument('--file', required=True, help='File to upload')
    parser.add_argument('--key', required=True, help='S3 object key')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.file):
        print(f"Error: File not found: {args.file}")
        return 1
    
    success = upload_to_s3(args.bucket, args.file, args.key)
    return 0 if success else 1


if __name__ == "__main__":
    exit(main())
