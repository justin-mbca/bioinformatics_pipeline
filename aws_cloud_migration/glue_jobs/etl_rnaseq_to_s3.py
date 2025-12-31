"""
AWS Glue ETL Job for RNA-Seq Data Processing

This script processes raw RNA-Seq count data from the AWS Glue Data Catalog,
filters low-quality reads, and writes the results to S3 in Parquet format.
"""

import sys
from typing import Dict, Any
from awsglue.transforms import *
from awsglue.utils import getResolvedOptions
from pyspark.context import SparkContext
from awsglue.context import GlueContext
from awsglue.job import Job
from pyspark.sql import DataFrame
from pyspark.sql.functions import col
import logging

# Configure logging
logger = logging.getLogger()
logger.setLevel(logging.INFO)


def initialize_glue_job(args: Dict[str, str]) -> tuple:
    """
    Initialize Glue job context and Spark session.
    
    Args:
        args: Dictionary containing job arguments
        
    Returns:
        Tuple of (SparkContext, GlueContext, Job)
    """
    sc = SparkContext()
    glueContext = GlueContext(sc)
    spark = glueContext.spark_session
    job = Job(glueContext)
    job.init(args['JOB_NAME'], args)
    
    logger.info(f"Initialized Glue job: {args['JOB_NAME']}")
    return sc, glueContext, job


def read_from_catalog(glueContext: GlueContext, database: str, table: str) -> DataFrame:
    """
    Read data from AWS Glue Data Catalog.
    
    Args:
        glueContext: GlueContext instance
        database: Glue database name
        table: Glue table name
        
    Returns:
        PySpark DataFrame with the data
    """
    logger.info(f"Reading from catalog - Database: {database}, Table: {table}")
    
    try:
        datasource = glueContext.create_dynamic_frame.from_catalog(
            database=database,
            table_name=table,
            transformation_ctx="datasource"
        )
        
        df = datasource.toDF()
        logger.info(f"Successfully read {df.count()} rows from catalog")
        return df
        
    except Exception as e:
        logger.error(f"Error reading from catalog: {str(e)}")
        raise


def filter_low_quality_reads(df: DataFrame, count_threshold: int = 10) -> DataFrame:
    """
    Filter out low-quality reads based on count threshold.
    
    Args:
        df: Input DataFrame with count data
        count_threshold: Minimum count value to retain
        
    Returns:
        Filtered DataFrame
    """
    logger.info(f"Filtering reads with count > {count_threshold}")
    
    try:
        # Get all count columns (excluding gene_id and other metadata)
        count_columns = [c for c in df.columns if c not in ['gene_id', 'gene_name', 'chromosome']]
        
        # Create filter condition: at least one sample must have count > threshold
        filter_condition = None
        for col_name in count_columns:
            condition = col(col_name) > count_threshold
            if filter_condition is None:
                filter_condition = condition
            else:
                filter_condition = filter_condition | condition
        
        filtered_df = df.filter(filter_condition)
        
        initial_count = df.count()
        filtered_count = filtered_df.count()
        logger.info(f"Filtered from {initial_count} to {filtered_count} genes "
                   f"({100 * filtered_count / initial_count:.2f}% retained)")
        
        return filtered_df
        
    except Exception as e:
        logger.error(f"Error filtering data: {str(e)}")
        raise


def write_to_s3(df: DataFrame, output_path: str, glueContext: GlueContext) -> None:
    """
    Write DataFrame to S3 in Parquet format.
    
    Args:
        df: DataFrame to write
        output_path: S3 path for output
        glueContext: GlueContext instance
    """
    logger.info(f"Writing data to S3: {output_path}")
    
    try:
        # Convert DataFrame to DynamicFrame
        from awsglue.dynamicframe import DynamicFrame
        
        dynamic_frame = DynamicFrame.fromDF(df, glueContext, "filtered_data")
        
        # Write using Glue
        glueContext.write_dynamic_frame.from_options(
            frame=dynamic_frame,
            connection_type="s3",
            connection_options={"path": output_path},
            format="parquet",
            transformation_ctx="datasink"
        )
        
        logger.info(f"Successfully wrote data to {output_path}")
        
    except Exception as e:
        logger.error(f"Error writing to S3: {str(e)}")
        raise


def main() -> None:
    """Main ETL workflow execution."""
    try:
        # Get job parameters
        args = getResolvedOptions(sys.argv, [
            'JOB_NAME',
            'database_name',
            'table_name',
            'output_path',
            'count_threshold'
        ])
        
        logger.info("Starting RNA-Seq ETL job")
        logger.info(f"Parameters: {args}")
        
        # Initialize Glue job
        sc, glueContext, job = initialize_glue_job(args)
        
        # Read data from catalog
        df = read_from_catalog(
            glueContext,
            args['database_name'],
            args['table_name']
        )
        
        # Filter low-quality reads
        count_threshold = int(args.get('count_threshold', '10'))
        filtered_df = filter_low_quality_reads(df, count_threshold)
        
        # Write results to S3
        write_to_s3(filtered_df, args['output_path'], glueContext)
        
        # Commit job
        job.commit()
        logger.info("RNA-Seq ETL job completed successfully")
        
    except Exception as e:
        logger.error(f"Job failed with error: {str(e)}")
        raise


if __name__ == "__main__":
    main()
