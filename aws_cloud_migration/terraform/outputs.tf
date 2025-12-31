# Terraform Outputs for Bioinformatics Pipeline

# S3 Outputs
output "s3_bucket_name" {
  description = "Name of the S3 data lake bucket"
  value       = aws_s3_bucket.data_lake.id
}

output "s3_bucket_arn" {
  description = "ARN of the S3 data lake bucket"
  value       = aws_s3_bucket.data_lake.arn
}

# Glue Outputs
output "glue_database_name" {
  description = "Name of the Glue database"
  value       = aws_glue_catalog_database.rnaseq_db.name
}

output "glue_job_name" {
  description = "Name of the Glue ETL job"
  value       = aws_glue_job.rnaseq_etl.name
}

# EMR Outputs
output "emr_cluster_id" {
  description = "ID of the EMR cluster"
  value       = aws_emr_cluster.bioinformatics.id
}

output "emr_master_public_dns" {
  description = "Public DNS of EMR master node"
  value       = aws_emr_cluster.bioinformatics.master_public_dns
}

# Lambda Outputs
output "trigger_pipeline_lambda_arn" {
  description = "ARN of trigger pipeline Lambda function"
  value       = aws_lambda_function.trigger_pipeline.arn
}

output "data_quality_check_lambda_arn" {
  description = "ARN of data quality check Lambda function"
  value       = aws_lambda_function.data_quality_check.arn
}

# SNS Outputs
output "sns_topic_arn" {
  description = "ARN of SNS topic for notifications"
  value       = aws_sns_topic.data_quality.arn
}

# IAM Outputs
output "glue_role_arn" {
  description = "ARN of Glue IAM role"
  value       = aws_iam_role.glue_role.arn
}

output "emr_service_role_arn" {
  description = "ARN of EMR service role"
  value       = aws_iam_role.emr_service_role.arn
}

output "emr_instance_profile_arn" {
  description = "ARN of EMR instance profile"
  value       = aws_iam_instance_profile.emr_instance_profile.arn
}

# Summary Output
output "deployment_summary" {
  description = "Summary of deployed resources"
  value = {
    region        = var.aws_region
    environment   = var.environment
    s3_bucket     = aws_s3_bucket.data_lake.id
    glue_database = aws_glue_catalog_database.rnaseq_db.name
    emr_cluster   = aws_emr_cluster.bioinformatics.id
    lambda_functions = [
      aws_lambda_function.trigger_pipeline.function_name,
      aws_lambda_function.data_quality_check.function_name
    ]
  }
}
