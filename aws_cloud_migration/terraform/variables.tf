# Terraform Variables for Bioinformatics Pipeline

variable "aws_region" {
  description = "AWS region for resources"
  type        = string
  default     = "us-east-1"
}

variable "environment" {
  description = "Environment name (dev, staging, production)"
  type        = string
  default     = "production"
}

variable "s3_bucket_name" {
  description = "Name of the S3 bucket for data lake"
  type        = string
  default     = "bioinformatics-pipeline"
}

variable "glue_database_name" {
  description = "Name of the Glue database"
  type        = string
  default     = "rnaseq_db"
}

variable "emr_instance_type" {
  description = "Instance type for EMR core nodes"
  type        = string
  default     = "r5.xlarge"
}

variable "emr_core_instance_count" {
  description = "Number of core instances for EMR cluster"
  type        = number
  default     = 3
}

variable "lambda_memory_size" {
  description = "Memory size (MB) for Lambda functions"
  type        = number
  default     = 512
}

variable "lambda_timeout" {
  description = "Timeout (seconds) for Lambda functions"
  type        = number
  default     = 300
}

variable "vpc_id" {
  description = "VPC ID for resources"
  type        = string
  default     = ""
}

variable "subnet_ids" {
  description = "List of subnet IDs"
  type        = list(string)
  default     = []
}

variable "sns_email" {
  description = "Email address for SNS notifications"
  type        = string
  default     = "data-team@example.com"
}

variable "enable_versioning" {
  description = "Enable S3 bucket versioning"
  type        = bool
  default     = true
}

variable "enable_logging" {
  description = "Enable CloudWatch logging"
  type        = bool
  default     = true
}

variable "tags" {
  description = "Additional tags for resources"
  type        = map(string)
  default     = {}
}
