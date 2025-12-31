"""
S3 Dataset Upload Script

This script uploads RNA-Seq datasets to S3 using boto3 with support for:
- Multi-part uploads for large files
- Progress tracking
- Metadata tagging
- Error handling and retries
"""

import boto3
import os
import sys
import logging
from pathlib import Path
from typing import Dict, Optional, List
from botocore.exceptions import ClientError
from botocore.config import Config
import hashlib
from datetime import datetime

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class S3DatasetUploader:
    """Class for uploading bioinformatics datasets to S3."""
    
    def __init__(
        self,
        bucket_name: str,
        region_name: str = 'us-east-1',
        multipart_threshold: int = 100 * 1024 * 1024  # 100 MB
    ):
        """
        Initialize S3 uploader.
        
        Args:
            bucket_name: S3 bucket name
            region_name: AWS region
            multipart_threshold: File size threshold for multipart upload (bytes)
        """
        self.bucket_name = bucket_name
        self.region_name = region_name
        self.multipart_threshold = multipart_threshold
        
        # Configure boto3 with retries
        config = Config(
            region_name=region_name,
            retries={'max_attempts': 3, 'mode': 'adaptive'}
        )
        
        self.s3_client = boto3.client('s3', config=config)
        self.s3_resource = boto3.resource('s3', config=config)
        
        logger.info(f"Initialized S3 uploader for bucket: {bucket_name}")
    
    def calculate_md5(self, file_path: str) -> str:
        """
        Calculate MD5 hash of file for integrity check.
        
        Args:
            file_path: Path to file
            
        Returns:
            MD5 hash as hex string
        """
        logger.info(f"Calculating MD5 hash for: {file_path}")
        
        hash_md5 = hashlib.md5()
        with open(file_path, "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                hash_md5.update(chunk)
        
        md5_hash = hash_md5.hexdigest()
        logger.info(f"MD5 hash: {md5_hash}")
        return md5_hash
    
    def upload_file(
        self,
        local_path: str,
        s3_key: str,
        metadata: Optional[Dict[str, str]] = None,
        tags: Optional[Dict[str, str]] = None,
        storage_class: str = 'STANDARD',
        show_progress: bool = True
    ) -> bool:
        """
        Upload file to S3 with metadata and tags.
        
        Args:
            local_path: Local file path
            s3_key: S3 object key (path in bucket)
            metadata: Custom metadata dictionary
            tags: Tags dictionary
            storage_class: S3 storage class
            show_progress: Whether to show progress bar
            
        Returns:
            True if successful, False otherwise
        """
        try:
            file_size = os.path.getsize(local_path)
            logger.info(f"Uploading {local_path} ({file_size:,} bytes) to s3://{self.bucket_name}/{s3_key}")
            
            # Prepare extra args
            extra_args = {
                'StorageClass': storage_class,
                'ServerSideEncryption': 'AES256'
            }
            
            # Add metadata
            if metadata is None:
                metadata = {}
            
            # Add default metadata
            metadata.update({
                'upload-date': datetime.now().isoformat(),
                'original-filename': os.path.basename(local_path),
                'file-size': str(file_size)
            })
            
            extra_args['Metadata'] = metadata
            
            # Add tags if provided
            if tags:
                tag_string = '&'.join([f"{k}={v}" for k, v in tags.items()])
                extra_args['Tagging'] = tag_string
            
            # Upload with progress callback
            if show_progress and file_size > self.multipart_threshold:
                class ProgressPercentage:
                    def __init__(self, filename, size):
                        self._filename = filename
                        self._size = size
                        self._seen_so_far = 0
                        
                    def __call__(self, bytes_amount):
                        self._seen_so_far += bytes_amount
                        percentage = (self._seen_so_far / self._size) * 100
                        sys.stdout.write(
                            f"\r{self._filename}: {self._seen_so_far:,} / {self._size:,} bytes "
                            f"({percentage:.2f}%)"
                        )
                        sys.stdout.flush()
                
                callback = ProgressPercentage(os.path.basename(local_path), file_size)
            else:
                callback = None
            
            # Perform upload
            self.s3_client.upload_file(
                local_path,
                self.bucket_name,
                s3_key,
                ExtraArgs=extra_args,
                Callback=callback
            )
            
            if callback:
                print()  # New line after progress
            
            logger.info(f"Successfully uploaded: s3://{self.bucket_name}/{s3_key}")
            return True
            
        except ClientError as e:
            logger.error(f"Error uploading file: {str(e)}")
            return False
        except Exception as e:
            logger.error(f"Unexpected error uploading file: {str(e)}")
            return False
    
    def upload_directory(
        self,
        local_dir: str,
        s3_prefix: str,
        file_patterns: Optional[List[str]] = None,
        metadata: Optional[Dict[str, str]] = None,
        tags: Optional[Dict[str, str]] = None
    ) -> Dict[str, int]:
        """
        Upload entire directory to S3.
        
        Args:
            local_dir: Local directory path
            s3_prefix: S3 key prefix
            file_patterns: List of file patterns to include (e.g., ['*.csv', '*.txt'])
            metadata: Metadata to apply to all files
            tags: Tags to apply to all files
            
        Returns:
            Dictionary with upload statistics
        """
        logger.info(f"Uploading directory: {local_dir} to s3://{self.bucket_name}/{s3_prefix}")
        
        stats = {'success': 0, 'failed': 0, 'skipped': 0}
        
        # Walk through directory
        for root, dirs, files in os.walk(local_dir):
            for filename in files:
                # Check file patterns
                if file_patterns:
                    if not any(Path(filename).match(pattern) for pattern in file_patterns):
                        stats['skipped'] += 1
                        continue
                
                # Construct paths
                local_path = os.path.join(root, filename)
                rel_path = os.path.relpath(local_path, local_dir)
                s3_key = os.path.join(s3_prefix, rel_path).replace('\\', '/')
                
                # Upload file
                success = self.upload_file(
                    local_path=local_path,
                    s3_key=s3_key,
                    metadata=metadata,
                    tags=tags
                )
                
                if success:
                    stats['success'] += 1
                else:
                    stats['failed'] += 1
        
        logger.info(f"Upload complete - Success: {stats['success']}, "
                   f"Failed: {stats['failed']}, Skipped: {stats['skipped']}")
        
        return stats
    
    def check_file_exists(self, s3_key: str) -> bool:
        """
        Check if file exists in S3.
        
        Args:
            s3_key: S3 object key
            
        Returns:
            True if exists, False otherwise
        """
        try:
            self.s3_client.head_object(Bucket=self.bucket_name, Key=s3_key)
            return True
        except ClientError:
            return False
    
    def list_files(self, prefix: str = '', max_keys: int = 1000) -> List[Dict[str, any]]:
        """
        List files in S3 bucket with prefix.
        
        Args:
            prefix: S3 key prefix to filter
            max_keys: Maximum number of keys to return
            
        Returns:
            List of file dictionaries with key, size, and last_modified
        """
        try:
            response = self.s3_client.list_objects_v2(
                Bucket=self.bucket_name,
                Prefix=prefix,
                MaxKeys=max_keys
            )
            
            if 'Contents' not in response:
                return []
            
            files = []
            for obj in response['Contents']:
                files.append({
                    'key': obj['Key'],
                    'size': obj['Size'],
                    'last_modified': obj['LastModified'],
                    'storage_class': obj.get('StorageClass', 'STANDARD')
                })
            
            logger.info(f"Found {len(files)} files with prefix: {prefix}")
            return files
            
        except ClientError as e:
            logger.error(f"Error listing files: {str(e)}")
            return []


def main():
    """Main upload workflow with example usage."""
    
    # Configuration
    BUCKET_NAME = os.environ.get('S3_BUCKET_NAME', 'bioinformatics-pipeline')
    REGION = os.environ.get('AWS_REGION', 'us-east-1')
    
    # Initialize uploader
    uploader = S3DatasetUploader(
        bucket_name=BUCKET_NAME,
        region_name=REGION
    )
    
    # Example 1: Upload single file
    logger.info("=== Example 1: Upload Single File ===")
    
    metadata = {
        'study': 'lung-cancer-2024',
        'data-type': 'rnaseq',
        'sample-count': '100',
        'processing-pipeline': 'v2.1.0'
    }
    
    tags = {
        'Project': 'bioinformatics-pipeline',
        'DataType': 'rnaseq',
        'Classification': 'sensitive',
        'Owner': 'data-science-team',
        'Environment': 'production'
    }
    
    # Uncomment to upload a real file
    # uploader.upload_file(
    #     local_path='path/to/local/counts.csv',
    #     s3_key='raw/counts/lung_cancer/cohort_a_counts.csv',
    #     metadata=metadata,
    #     tags=tags,
    #     storage_class='STANDARD'
    # )
    
    # Example 2: Upload directory
    logger.info("=== Example 2: Upload Directory ===")
    
    # Uncomment to upload a directory
    # stats = uploader.upload_directory(
    #     local_dir='path/to/local/data',
    #     s3_prefix='raw/counts/lung_cancer/',
    #     file_patterns=['*.csv', '*.tsv'],
    #     metadata=metadata,
    #     tags=tags
    # )
    
    # Example 3: List uploaded files
    logger.info("=== Example 3: List Uploaded Files ===")
    
    files = uploader.list_files(prefix='raw/counts/', max_keys=10)
    
    for file_info in files:
        logger.info(f"File: {file_info['key']}, "
                   f"Size: {file_info['size']:,} bytes, "
                   f"Modified: {file_info['last_modified']}")
    
    logger.info("Upload script completed")


if __name__ == "__main__":
    main()
