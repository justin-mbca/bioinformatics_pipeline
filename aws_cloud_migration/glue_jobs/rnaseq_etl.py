"""
AWS Glue ETL Job - Process RNA-Seq Data
"""
import sys
from awsglue.transforms import *
from awsglue.utils import getResolvedOptions
from pyspark.context import SparkContext
from awsglue.context import GlueContext
from awsglue.job import Job


def main():
    """Main Glue job function"""
    args = getResolvedOptions(sys.argv, ['JOB_NAME', 'S3_INPUT_PATH', 'S3_OUTPUT_PATH'])
    
    sc = SparkContext()
    glueContext = GlueContext(sc)
    spark = glueContext.spark_session
    job = Job(glueContext)
    job.init(args['JOB_NAME'], args)
    
    # Read data from S3
    input_path = args['S3_INPUT_PATH']
    output_path = args['S3_OUTPUT_PATH']
    
    print(f"Reading data from: {input_path}")
    df = spark.read.csv(input_path, header=True, inferSchema=True)
    
    # Process data
    print(f"Processing {df.count()} rows")
    
    # Write processed data
    print(f"Writing data to: {output_path}")
    df.write.mode('overwrite').parquet(output_path)
    
    job.commit()
    print("Glue job completed successfully")


if __name__ == "__main__":
    main()
