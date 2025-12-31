# Main Terraform Configuration for Bioinformatics Pipeline
# This file configures the AWS provider and backend for state management

terraform {
  required_version = ">= 1.0"
  
  required_providers {
    aws = {
      source  = "hashicorp/aws"
      version = "~> 5.0"
    }
  }
  
  # S3 backend for state management
  backend "s3" {
    bucket         = "bioinformatics-terraform-state"
    key            = "pipeline/terraform.tfstate"
    region         = "us-east-1"
    encrypt        = true
    dynamodb_table = "terraform-state-lock"
  }
}

# AWS Provider Configuration
provider "aws" {
  region = var.aws_region
  
  default_tags {
    tags = {
      Project     = "bioinformatics-pipeline"
      Environment = var.environment
      ManagedBy   = "Terraform"
      Owner       = "data-science-team"
    }
  }
}

# Data sources
data "aws_caller_identity" "current" {}
data "aws_region" "current" {}
