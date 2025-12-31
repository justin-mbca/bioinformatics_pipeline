# Terraform configuration for AWS bioinformatics infrastructure

terraform {
  required_version = ">= 1.0"

  required_providers {
    aws = {
      source  = "hashicorp/aws"
      version = "~> 5.0"
    }
  }
}

provider "aws" {
  region = var.aws_region
}

variable "aws_region" {
  description = "AWS region"
  type        = string
  default     = "us-east-1"
}

variable "project_name" {
  description = "Project name"
  type        = string
  default     = "bioinformatics"
}

# S3 bucket for data lake
resource "aws_s3_bucket" "data_lake" {
  bucket = "${var.project_name}-data-lake"

  tags = {
    Name    = "Bioinformatics Data Lake"
    Project = var.project_name
  }
}

# S3 bucket versioning
resource "aws_s3_bucket_versioning" "data_lake" {
  bucket = aws_s3_bucket.data_lake.id

  versioning_configuration {
    status = "Enabled"
  }
}

# Output values
output "data_lake_bucket" {
  value       = aws_s3_bucket.data_lake.id
  description = "S3 data lake bucket name"
}
