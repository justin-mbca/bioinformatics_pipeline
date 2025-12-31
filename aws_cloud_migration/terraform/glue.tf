# AWS Glue Resources

# Glue Database
resource "aws_glue_catalog_database" "rnaseq_db" {
  name        = var.glue_database_name
  description = "RNA-Seq data catalog database"

  location_uri = "s3://${aws_s3_bucket.data_lake.id}/catalog/${var.glue_database_name}/"
}

# IAM Role for Glue
resource "aws_iam_role" "glue_role" {
  name = "bioinformatics-glue-role"

  assume_role_policy = jsonencode({
    Version = "2012-10-17"
    Statement = [
      {
        Action = "sts:AssumeRole"
        Effect = "Allow"
        Principal = {
          Service = "glue.amazonaws.com"
        }
      }
    ]
  })
}

# Glue Service Policy
resource "aws_iam_role_policy_attachment" "glue_service" {
  role       = aws_iam_role.glue_role.name
  policy_arn = "arn:aws:iam::aws:policy/service-role/AWSGlueServiceRole"
}

# S3 Access Policy for Glue
resource "aws_iam_role_policy" "glue_s3_access" {
  name = "glue-s3-access"
  role = aws_iam_role.glue_role.id

  policy = jsonencode({
    Version = "2012-10-17"
    Statement = [
      {
        Effect = "Allow"
        Action = [
          "s3:GetObject",
          "s3:PutObject",
          "s3:DeleteObject",
          "s3:ListBucket"
        ]
        Resource = [
          aws_s3_bucket.data_lake.arn,
          "${aws_s3_bucket.data_lake.arn}/*"
        ]
      }
    ]
  })
}

# Glue ETL Job
resource "aws_glue_job" "rnaseq_etl" {
  name     = "rnaseq-etl-pipeline"
  role_arn = aws_iam_role.glue_role.arn

  command {
    name            = "glueetl"
    script_location = "s3://${aws_s3_bucket.data_lake.id}/glue-scripts/etl_rnaseq_to_s3.py"
    python_version  = "3"
  }

  default_arguments = {
    "--job-language"                     = "python"
    "--job-bookmark-option"              = "job-bookmark-enable"
    "--enable-metrics"                   = "true"
    "--enable-continuous-cloudwatch-log" = "true"
    "--enable-spark-ui"                  = "true"
    "--spark-event-logs-path"            = "s3://${aws_s3_bucket.data_lake.id}/glue-logs/spark-ui/"
    "--TempDir"                          = "s3://${aws_s3_bucket.data_lake.id}/glue-temp/"
    "--database_name"                    = var.glue_database_name
    "--table_name"                       = "raw_counts"
    "--output_path"                      = "s3://${aws_s3_bucket.data_lake.id}/processed/normalized/"
    "--count_threshold"                  = "10"
  }

  max_retries  = 1
  timeout      = 2880
  glue_version = "4.0"
  worker_type  = "G.1X"
  number_of_workers = 5

  execution_property {
    max_concurrent_runs = 3
  }
}

# Glue Crawler
resource "aws_glue_crawler" "raw_counts" {
  database_name = aws_glue_catalog_database.rnaseq_db.name
  name          = "raw-counts-crawler"
  role          = aws_iam_role.glue_role.arn

  s3_target {
    path = "s3://${aws_s3_bucket.data_lake.id}/raw/counts/"
  }

  schedule = "cron(0 2 * * ? *)"  # Daily at 2 AM

  schema_change_policy {
    delete_behavior = "LOG"
    update_behavior = "UPDATE_IN_DATABASE"
  }
}
