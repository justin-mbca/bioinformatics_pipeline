# AWS EMR Cluster Resources

# EMR Service Role
resource "aws_iam_role" "emr_service_role" {
  name = "bioinformatics-emr-service-role"

  assume_role_policy = jsonencode({
    Version = "2012-10-17"
    Statement = [
      {
        Action = "sts:AssumeRole"
        Effect = "Allow"
        Principal = {
          Service = "elasticmapreduce.amazonaws.com"
        }
      }
    ]
  })
}

resource "aws_iam_role_policy_attachment" "emr_service_policy" {
  role       = aws_iam_role.emr_service_role.name
  policy_arn = "arn:aws:iam::aws:policy/service-role/AmazonElasticMapReduceRole"
}

# EMR EC2 Instance Profile Role
resource "aws_iam_role" "emr_ec2_role" {
  name = "bioinformatics-emr-ec2-role"

  assume_role_policy = jsonencode({
    Version = "2012-10-17"
    Statement = [
      {
        Action = "sts:AssumeRole"
        Effect = "Allow"
        Principal = {
          Service = "ec2.amazonaws.com"
        }
      }
    ]
  })
}

resource "aws_iam_role_policy_attachment" "emr_ec2_policy" {
  role       = aws_iam_role.emr_ec2_role.name
  policy_arn = "arn:aws:iam::aws:policy/service-role/AmazonElasticMapReduceforEC2Role"
}

# S3 Access Policy for EMR
resource "aws_iam_role_policy" "emr_s3_access" {
  name = "emr-s3-access"
  role = aws_iam_role.emr_ec2_role.id

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

resource "aws_iam_instance_profile" "emr_instance_profile" {
  name = "bioinformatics-emr-instance-profile"
  role = aws_iam_role.emr_ec2_role.name
}

# EMR Security Configuration
resource "aws_emr_security_configuration" "bioinformatics" {
  name = "bioinformatics-emr-security-config"

  configuration = jsonencode({
    EncryptionConfiguration = {
      EnableInTransitEncryption = false
      EnableAtRestEncryption    = true
      AtRestEncryptionConfiguration = {
        S3EncryptionConfiguration = {
          EncryptionMode = "SSE-S3"
        }
        LocalDiskEncryptionConfiguration = {
          EncryptionKeyProviderType = "AwsKms"
          AwsKmsKey                 = "alias/aws/elasticmapreduce"
        }
      }
    }
  })
}

# EMR Cluster
resource "aws_emr_cluster" "bioinformatics" {
  name          = "bioinformatics-deseq2-cluster"
  release_label = "emr-6.10.0"
  applications  = ["Spark", "Hadoop", "Hive", "Livy"]

  service_role = aws_iam_role.emr_service_role.arn

  ec2_attributes {
    instance_profile                  = aws_iam_instance_profile.emr_instance_profile.arn
    emr_managed_master_security_group = aws_security_group.emr_master.id
    emr_managed_slave_security_group  = aws_security_group.emr_core.id
  }

  master_instance_group {
    instance_type  = "m5.xlarge"
    instance_count = 1

    ebs_config {
      size                 = 100
      type                 = "gp3"
      volumes_per_instance = 1
    }
  }

  core_instance_group {
    instance_type  = var.emr_instance_type
    instance_count = var.emr_core_instance_count

    ebs_config {
      size                 = 200
      type                 = "gp3"
      volumes_per_instance = 1
    }

    autoscaling_policy = jsonencode({
      Constraints = {
        MinCapacity = 2
        MaxCapacity = 10
      }
      Rules = [
        {
          Name = "ScaleUpOnYARNMemory"
          Action = {
            SimpleScalingPolicyConfiguration = {
              AdjustmentType    = "CHANGE_IN_CAPACITY"
              ScalingAdjustment = 2
              CoolDown          = 300
            }
          }
          Trigger = {
            CloudWatchAlarmDefinition = {
              ComparisonOperator = "GREATER_THAN"
              EvaluationPeriods  = 2
              MetricName         = "YARNMemoryAvailablePercentage"
              Namespace          = "AWS/ElasticMapReduce"
              Period             = 300
              Statistic          = "AVERAGE"
              Threshold          = 75.0
              Unit               = "PERCENT"
            }
          }
        }
      ]
    })
  }

  bootstrap_action {
    path = "s3://${aws_s3_bucket.data_lake.id}/bootstrap/bootstrap_actions.sh"
    name = "Install Bioinformatics Packages"
  }

  configurations_json = jsonencode([
    {
      Classification = "spark"
      Properties = {
        maximizeResourceAllocation = "true"
      }
    },
    {
      Classification = "spark-defaults"
      Properties = {
        "spark.dynamicAllocation.enabled" = "true"
        "spark.shuffle.service.enabled"   = "true"
        "spark.sql.adaptive.enabled"      = "true"
        "spark.serializer"                = "org.apache.spark.serializer.KryoSerializer"
      }
    }
  ])

  log_uri = "s3://${aws_s3_bucket.data_lake.id}/logs/emr/"

  keep_job_flow_alive_when_no_steps = true
  termination_protection            = false

  auto_termination_policy {
    idle_timeout = 3600
  }

  security_configuration = aws_emr_security_configuration.bioinformatics.name

  tags = merge(
    var.tags,
    {
      Name = "bioinformatics-deseq2-cluster"
    }
  )
}

# Security Groups
resource "aws_security_group" "emr_master" {
  name        = "emr-master-sg"
  description = "Security group for EMR master node"
  vpc_id      = var.vpc_id != "" ? var.vpc_id : data.aws_vpc.default.id

  egress {
    from_port   = 0
    to_port     = 0
    protocol    = "-1"
    cidr_blocks = ["0.0.0.0/0"]
  }

  tags = {
    Name = "emr-master-sg"
  }
}

resource "aws_security_group" "emr_core" {
  name        = "emr-core-sg"
  description = "Security group for EMR core nodes"
  vpc_id      = var.vpc_id != "" ? var.vpc_id : data.aws_vpc.default.id

  egress {
    from_port   = 0
    to_port     = 0
    protocol    = "-1"
    cidr_blocks = ["0.0.0.0/0"]
  }

  tags = {
    Name = "emr-core-sg"
  }
}

# Data source for default VPC
data "aws_vpc" "default" {
  default = true
}
