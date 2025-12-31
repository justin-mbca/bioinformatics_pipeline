"""
AWS Glue Data Catalog Setup Script

This script creates the Glue Data Catalog database and table schema
for RNA-Seq count data using boto3.
"""

import boto3
import json
import logging
from typing import Dict, List, Any, Optional
from botocore.exceptions import ClientError

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class GlueCatalogManager:
    """Manager class for AWS Glue Data Catalog operations."""
    
    def __init__(self, region_name: str = 'us-east-1'):
        """
        Initialize Glue catalog manager.
        
        Args:
            region_name: AWS region name
        """
        self.glue_client = boto3.client('glue', region_name=region_name)
        self.region_name = region_name
        logger.info(f"Initialized Glue catalog manager for region: {region_name}")
    
    def create_database(self, database_name: str, description: str = '') -> bool:
        """
        Create a Glue database.
        
        Args:
            database_name: Name of the database to create
            description: Optional description for the database
            
        Returns:
            True if successful, False otherwise
        """
        try:
            response = self.glue_client.create_database(
                DatabaseInput={
                    'Name': database_name,
                    'Description': description or f'RNA-Seq data catalog database',
                    'LocationUri': f's3://bioinformatics-pipeline/catalog/{database_name}/',
                    'Parameters': {
                        'classification': 'genomics',
                        'data_type': 'rnaseq',
                        'created_by': 'glue_catalog_setup_script'
                    }
                }
            )
            logger.info(f"Successfully created database: {database_name}")
            return True
            
        except ClientError as e:
            if e.response['Error']['Code'] == 'AlreadyExistsException':
                logger.warning(f"Database {database_name} already exists")
                return True
            else:
                logger.error(f"Error creating database: {str(e)}")
                return False
    
    def get_rnaseq_table_schema(self) -> List[Dict[str, str]]:
        """
        Define the table schema for raw RNA-Seq counts.
        
        Returns:
            List of column definitions
        """
        return [
            {'Name': 'gene_id', 'Type': 'string', 'Comment': 'Ensembl gene identifier'},
            {'Name': 'gene_name', 'Type': 'string', 'Comment': 'HGNC gene symbol'},
            {'Name': 'chromosome', 'Type': 'string', 'Comment': 'Chromosome location'},
            {'Name': 'start_position', 'Type': 'int', 'Comment': 'Gene start position'},
            {'Name': 'end_position', 'Type': 'int', 'Comment': 'Gene end position'},
            {'Name': 'strand', 'Type': 'string', 'Comment': 'Strand orientation (+/-)'},
            {'Name': 'gene_length', 'Type': 'int', 'Comment': 'Gene length in bases'},
            {'Name': 'sample_01', 'Type': 'int', 'Comment': 'Read count for sample 01'},
            {'Name': 'sample_02', 'Type': 'int', 'Comment': 'Read count for sample 02'},
            {'Name': 'sample_03', 'Type': 'int', 'Comment': 'Read count for sample 03'},
            {'Name': 'sample_04', 'Type': 'int', 'Comment': 'Read count for sample 04'},
            {'Name': 'sample_05', 'Type': 'int', 'Comment': 'Read count for sample 05'},
            {'Name': 'sample_06', 'Type': 'int', 'Comment': 'Read count for sample 06'},
        ]
    
    def create_table(
        self,
        database_name: str,
        table_name: str,
        s3_location: str,
        schema: Optional[List[Dict[str, str]]] = None,
        input_format: str = 'org.apache.hadoop.mapred.TextInputFormat',
        output_format: str = 'org.apache.hadoop.hive.ql.io.HiveIgnoreKeyTextOutputFormat',
        serde_info: Optional[Dict[str, Any]] = None
    ) -> bool:
        """
        Create a Glue table with specified schema.
        
        Args:
            database_name: Name of the database
            table_name: Name of the table to create
            s3_location: S3 location of the data
            schema: Column definitions
            input_format: Hadoop input format class
            output_format: Hadoop output format class
            serde_info: SerDe information
            
        Returns:
            True if successful, False otherwise
        """
        if schema is None:
            schema = self.get_rnaseq_table_schema()
        
        if serde_info is None:
            serde_info = {
                'SerializationLibrary': 'org.apache.hadoop.hive.serde2.lazy.LazySimpleSerDe',
                'Parameters': {
                    'field.delim': ',',
                    'skip.header.line.count': '1'
                }
            }
        
        table_input = {
            'Name': table_name,
            'Description': 'Raw RNA-Seq count data',
            'StorageDescriptor': {
                'Columns': schema,
                'Location': s3_location,
                'InputFormat': input_format,
                'OutputFormat': output_format,
                'Compressed': False,
                'SerdeInfo': serde_info,
                'StoredAsSubDirectories': False
            },
            'PartitionKeys': [],
            'TableType': 'EXTERNAL_TABLE',
            'Parameters': {
                'classification': 'csv',
                'compressionType': 'none',
                'typeOfData': 'file'
            }
        }
        
        try:
            response = self.glue_client.create_table(
                DatabaseName=database_name,
                TableInput=table_input
            )
            logger.info(f"Successfully created table: {database_name}.{table_name}")
            return True
            
        except ClientError as e:
            if e.response['Error']['Code'] == 'AlreadyExistsException':
                logger.warning(f"Table {database_name}.{table_name} already exists")
                return True
            else:
                logger.error(f"Error creating table: {str(e)}")
                return False
    
    def create_parquet_table(
        self,
        database_name: str,
        table_name: str,
        s3_location: str,
        schema: Optional[List[Dict[str, str]]] = None
    ) -> bool:
        """
        Create a Glue table for Parquet format data.
        
        Args:
            database_name: Name of the database
            table_name: Name of the table to create
            s3_location: S3 location of the Parquet data
            schema: Column definitions
            
        Returns:
            True if successful, False otherwise
        """
        if schema is None:
            schema = self.get_rnaseq_table_schema()
        
        serde_info = {
            'SerializationLibrary': 'org.apache.hadoop.hive.ql.io.parquet.serde.ParquetHiveSerDe',
            'Parameters': {
                'serialization.format': '1'
            }
        }
        
        return self.create_table(
            database_name=database_name,
            table_name=table_name,
            s3_location=s3_location,
            schema=schema,
            input_format='org.apache.hadoop.hive.ql.io.parquet.MapredParquetInputFormat',
            output_format='org.apache.hadoop.hive.ql.io.parquet.MapredParquetOutputFormat',
            serde_info=serde_info
        )
    
    def list_databases(self) -> List[str]:
        """
        List all Glue databases.
        
        Returns:
            List of database names
        """
        try:
            response = self.glue_client.get_databases()
            databases = [db['Name'] for db in response.get('DatabaseList', [])]
            logger.info(f"Found {len(databases)} databases")
            return databases
            
        except ClientError as e:
            logger.error(f"Error listing databases: {str(e)}")
            return []
    
    def list_tables(self, database_name: str) -> List[str]:
        """
        List all tables in a database.
        
        Args:
            database_name: Name of the database
            
        Returns:
            List of table names
        """
        try:
            response = self.glue_client.get_tables(DatabaseName=database_name)
            tables = [table['Name'] for table in response.get('TableList', [])]
            logger.info(f"Found {len(tables)} tables in database {database_name}")
            return tables
            
        except ClientError as e:
            logger.error(f"Error listing tables: {str(e)}")
            return []


def main():
    """Main setup workflow."""
    # Configuration
    DATABASE_NAME = 'rnaseq_db'
    RAW_TABLE_NAME = 'raw_counts'
    PROCESSED_TABLE_NAME = 'processed_counts'
    RAW_S3_LOCATION = 's3://bioinformatics-pipeline/raw/counts/'
    PROCESSED_S3_LOCATION = 's3://bioinformatics-pipeline/processed/normalized/'
    REGION = 'us-east-1'
    
    logger.info("Starting Glue Data Catalog setup")
    
    # Initialize manager
    manager = GlueCatalogManager(region_name=REGION)
    
    # Create database
    success = manager.create_database(
        database_name=DATABASE_NAME,
        description='RNA-Seq analysis data catalog'
    )
    
    if not success:
        logger.error("Failed to create database. Exiting.")
        return
    
    # Create raw counts table (CSV format)
    success = manager.create_table(
        database_name=DATABASE_NAME,
        table_name=RAW_TABLE_NAME,
        s3_location=RAW_S3_LOCATION
    )
    
    if not success:
        logger.error("Failed to create raw counts table")
    
    # Create processed counts table (Parquet format)
    success = manager.create_parquet_table(
        database_name=DATABASE_NAME,
        table_name=PROCESSED_TABLE_NAME,
        s3_location=PROCESSED_S3_LOCATION
    )
    
    if not success:
        logger.error("Failed to create processed counts table")
    
    # List created resources
    logger.info("\n=== Created Resources ===")
    databases = manager.list_databases()
    logger.info(f"Databases: {databases}")
    
    if DATABASE_NAME in databases:
        tables = manager.list_tables(DATABASE_NAME)
        logger.info(f"Tables in {DATABASE_NAME}: {tables}")
    
    logger.info("Glue Data Catalog setup completed")


if __name__ == "__main__":
    main()
