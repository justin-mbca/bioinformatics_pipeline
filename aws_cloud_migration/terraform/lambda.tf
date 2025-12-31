# AWS Lambda Functions

# SNS Topic for notifications
resource "aws_sns_topic" "data_quality" {
  name = "bioinformatics-data-quality"
  display_name = "Bioinformatics Data Quality Notifications"
}

resource "aws_sns_topic_subscription" "data_quality_email" {
  topic_arn = aws_sns_topic.data_quality.arn
  protocol  = "email"
  endpoint  = var.sns_email
}

# IAM Role for Lambda
resource "aws_iam_role" "lambda_role" {
  name = "bioinformatics-lambda-role"

  assume_role_policy = jsonencode({
    Version = "2012-10-17"
    Statement = [
      {
        Action = "sts:AssumeRole"
        Effect = "Allow"
        Principal = {
          Service = "lambda.amazonaws.com"
        }
      }
    ]
  })
}

# Lambda Basic Execution Policy
resource "aws_iam_role_policy_attachment" "lambda_basic" {
  role       = aws_iam_role.lambda_role.name
  policy_arn = "arn:aws:iam::aws:policy/service-role/AWSLambdaBasicExecutionRole"
}

# S3 Access Policy for Lambda
resource "aws_iam_role_policy" "lambda_s3_access" {
  name = "lambda-s3-access"
  role = aws_iam_role.lambda_role.id

  policy = jsonencode({
    Version = "2012-10-17"
    Statement = [
      {
        Effect = "Allow"
        Action = [
          "s3:GetObject",
          "s3:PutObject",
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

# Glue Access Policy for Lambda
resource "aws_iam_role_policy" "lambda_glue_access" {
  name = "lambda-glue-access"
  role = aws_iam_role.lambda_role.id

  policy = jsonencode({
    Version = "2012-10-17"
    Statement = [
      {
        Effect = "Allow"
        Action = [
          "glue:StartJobRun",
          "glue:GetJobRun",
          "glue:GetJobRuns",
          "glue:BatchStopJobRun"
        ]
        Resource = "*"
      }
    ]
  })
}

# SNS Publish Policy for Lambda
resource "aws_iam_role_policy" "lambda_sns_publish" {
  name = "lambda-sns-publish"
  role = aws_iam_role.lambda_role.id

  policy = jsonencode({
    Version = "2012-10-17"
    Statement = [
      {
        Effect = "Allow"
        Action = [
          "sns:Publish"
        ]
        Resource = aws_sns_topic.data_quality.arn
      }
    ]
  })
}

# Lambda Function: Trigger Pipeline
resource "aws_lambda_function" "trigger_pipeline" {
  filename      = "${path.module}/../lambda_functions/trigger_pipeline.zip"
  function_name = "bioinformatics-trigger-pipeline"
  role          = aws_iam_role.lambda_role.arn
  handler       = "trigger_pipeline.lambda_handler"
  runtime       = "python3.9"
  timeout       = var.lambda_timeout
  memory_size   = var.lambda_memory_size

  environment {
    variables = {
      GLUE_JOB_NAME   = aws_glue_job.rnaseq_etl.name
      GLUE_DATABASE   = var.glue_database_name
      GLUE_TABLE      = "raw_counts"
      OUTPUT_BUCKET   = var.s3_bucket_name
      COUNT_THRESHOLD = "10"
    }
  }

  lifecycle {
    ignore_changes = [filename]
  }
}

# Lambda Function: Data Quality Check
resource "aws_lambda_function" "data_quality_check" {
  filename      = "${path.module}/../lambda_functions/data_quality_check.zip"
  function_name = "bioinformatics-data-quality-check"
  role          = aws_iam_role.lambda_role.arn
  handler       = "data_quality_check.lambda_handler"
  runtime       = "python3.9"
  timeout       = var.lambda_timeout
  memory_size   = var.lambda_memory_size

  environment {
    variables = {
      SNS_TOPIC_ARN       = aws_sns_topic.data_quality.arn
      MAX_MISSING_PERCENT = "10.0"
      MIN_SAMPLE_COUNT    = "3"
      REQUIRED_COLUMNS    = "gene_id,gene_name"
    }
  }

  lifecycle {
    ignore_changes = [filename]
  }
}

# S3 Event Notification for Trigger Pipeline Lambda
resource "aws_s3_bucket_notification" "trigger_pipeline" {
  bucket = aws_s3_bucket.data_lake.id

  lambda_function {
    lambda_function_arn = aws_lambda_function.trigger_pipeline.arn
    events              = ["s3:ObjectCreated:*"]
    filter_prefix       = "raw/counts/"
    filter_suffix       = ".csv"
  }

  depends_on = [aws_lambda_permission.allow_s3_trigger]
}

# Lambda Permission for S3 to invoke function
resource "aws_lambda_permission" "allow_s3_trigger" {
  statement_id  = "AllowExecutionFromS3"
  action        = "lambda:InvokeFunction"
  function_name = aws_lambda_function.trigger_pipeline.function_name
  principal     = "s3.amazonaws.com"
  source_arn    = aws_s3_bucket.data_lake.arn
}

# CloudWatch Log Groups
resource "aws_cloudwatch_log_group" "trigger_pipeline" {
  name              = "/aws/lambda/${aws_lambda_function.trigger_pipeline.function_name}"
  retention_in_days = 30
}

resource "aws_cloudwatch_log_group" "data_quality_check" {
  name              = "/aws/lambda/${aws_lambda_function.data_quality_check.function_name}"
  retention_in_days = 30
}
