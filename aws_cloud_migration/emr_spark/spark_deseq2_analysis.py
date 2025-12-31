"""
EMR Spark Job for Differential Expression Analysis

This PySpark script performs DESeq2-style differential expression analysis
on RNA-Seq count data stored in S3. It implements normalization, fold change
calculation, and statistical testing at scale.
"""

from pyspark.sql import SparkSession, DataFrame
from pyspark.sql.functions import (
    col, log, exp, sum as _sum, mean as _mean, stddev, count, lit,
    when, abs as _abs, pow as _pow, sqrt, array, stddev_samp, row_number,
    count as _count, min as _min
)
from pyspark.sql.window import Window
from pyspark.sql.types import DoubleType, StringType, StructType, StructField, IntegerType
from typing import List, Tuple
import logging
import sys

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class DESeq2SparkAnalyzer:
    """
    Spark-based implementation of DESeq2-style differential expression analysis.
    
    Implements key DESeq2 methods:
    - Size factor estimation (median of ratios)
    - Normalization
    - Log2 fold change calculation
    - Statistical significance testing
    """
    
    def __init__(self, spark: SparkSession):
        """
        Initialize analyzer with Spark session.
        
        Args:
            spark: Active SparkSession
        """
        self.spark = spark
        logger.info("Initialized DESeq2SparkAnalyzer")
    
    def read_count_data(self, input_path: str) -> DataFrame:
        """
        Read count data from S3 Parquet files.
        
        Args:
            input_path: S3 path to Parquet files
            
        Returns:
            DataFrame with count data
        """
        logger.info(f"Reading count data from: {input_path}")
        df = self.spark.read.parquet(input_path)
        logger.info(f"Loaded {df.count()} genes")
        return df
    
    def calculate_size_factors(
        self,
        df: DataFrame,
        sample_columns: List[str]
    ) -> dict:
        """
        Calculate DESeq2-style size factors using median of ratios method.
        
        Args:
            df: DataFrame with count data
            sample_columns: List of sample column names
            
        Returns:
            Dictionary mapping sample names to size factors
        """
        logger.info("Calculating size factors (median of ratios method)")
        
        # Step 1: Calculate geometric mean for each gene across samples
        # Log transform and average, then exp to get geometric mean
        log_sum_expr = " + ".join([f"log(NULLIF({col}, 0))" for col in sample_columns])
        geometric_mean_expr = f"exp(({log_sum_expr}) / {len(sample_columns)})"
        
        df_with_geo_mean = df.selectExpr(
            "*",
            f"{geometric_mean_expr} as geo_mean"
        ).filter("geo_mean > 0")
        
        # Step 2: Calculate ratio of each count to geometric mean
        size_factors = {}
        
        for sample_col in sample_columns:
            # Calculate ratios for this sample
            df_ratios = df_with_geo_mean.selectExpr(
                f"{sample_col} / geo_mean as ratio"
            ).filter("ratio > 0 AND ratio IS NOT NULL")
            
            # Get median ratio (approximate using percentile_approx for large datasets)
            median_ratio = df_ratios.stat.approxQuantile("ratio", [0.5], 0.01)[0]
            size_factors[sample_col] = median_ratio
            
            logger.info(f"Size factor for {sample_col}: {median_ratio:.4f}")
        
        return size_factors
    
    def normalize_counts(
        self,
        df: DataFrame,
        sample_columns: List[str],
        size_factors: dict
    ) -> DataFrame:
        """
        Normalize counts by dividing by size factors.
        
        Args:
            df: DataFrame with raw counts
            sample_columns: List of sample column names
            size_factors: Dictionary of size factors per sample
            
        Returns:
            DataFrame with normalized counts
        """
        logger.info("Normalizing counts by size factors")
        
        # Create expressions for normalized columns
        for sample_col in sample_columns:
            sf = size_factors[sample_col]
            df = df.withColumn(
                f"{sample_col}_norm",
                col(sample_col) / lit(sf)
            )
        
        return df
    
    def calculate_log2_fold_change(
        self,
        df: DataFrame,
        condition1_samples: List[str],
        condition2_samples: List[str],
        pseudocount: float = 1.0
    ) -> DataFrame:
        """
        Calculate log2 fold change between two conditions.
        
        Args:
            df: DataFrame with normalized counts
            condition1_samples: Sample columns for condition 1 (e.g., treatment)
            condition2_samples: Sample columns for condition 2 (e.g., control)
            pseudocount: Pseudocount to add before log transformation
            
        Returns:
            DataFrame with log2FC, baseMean, and group means
        """
        logger.info("Calculating log2 fold changes")
        
        # Calculate mean for each condition
        cond1_cols = [f"{s}_norm" for s in condition1_samples]
        cond2_cols = [f"{s}_norm" for s in condition2_samples]
        
        # Mean expression for condition 1
        cond1_mean_expr = " + ".join(cond1_cols) + f" / {len(cond1_cols)}"
        
        # Mean expression for condition 2
        cond2_mean_expr = " + ".join(cond2_cols) + f" / {len(cond2_cols)}"
        
        # Overall base mean
        all_norm_cols = cond1_cols + cond2_cols
        base_mean_expr = " + ".join(all_norm_cols) + f" / {len(all_norm_cols)}"
        
        df = df.selectExpr(
            "*",
            f"{cond1_mean_expr} as cond1_mean",
            f"{cond2_mean_expr} as cond2_mean",
            f"{base_mean_expr} as baseMean"
        )
        
        # Calculate log2 fold change
        df = df.withColumn(
            "log2FoldChange",
            log((col("cond1_mean") + lit(pseudocount)) / 
                (col("cond2_mean") + lit(pseudocount))) / log(lit(2.0))
        )
        
        return df
    
    def calculate_statistics(
        self,
        df: DataFrame,
        condition1_samples: List[str],
        condition2_samples: List[str]
    ) -> DataFrame:
        """
        Calculate statistical significance using Wald test approximation.
        
        Args:
            df: DataFrame with normalized counts and fold changes
            condition1_samples: Sample columns for condition 1
            condition2_samples: Sample columns for condition 2
            
        Returns:
            DataFrame with pvalue and padj columns
        """
        logger.info("Calculating statistical significance")
        
        # Get normalized column names
        cond1_cols = [f"{s}_norm" for s in condition1_samples]
        cond2_cols = [f"{s}_norm" for s in condition2_samples]
        
        # Calculate pooled standard error
        # This is a simplified version - production would use more sophisticated methods
        
        # Create arrays of values for each condition
        df = df.withColumn("cond1_values", array(*[col(c) for c in cond1_cols]))
        df = df.withColumn("cond2_values", array(*[col(c) for c in cond2_cols]))
        
        # Simplified variance calculation (would need UDF for proper implementation)
        # For now, use approximate method based on base mean and fold change
        df = df.withColumn(
            "wald_statistic",
            col("log2FoldChange") / (lit(1.0) / sqrt(col("baseMean") + lit(1.0)))
        )
        
        # Convert to p-value (approximate using standard normal)
        # In production, would use proper Wald test with negative binomial
        from pyspark.sql.functions import when, abs as _abs
        
        df = df.withColumn(
            "pvalue",
            when(col("baseMean") < 1, lit(1.0))
            .otherwise(
                # Approximate p-value from z-score
                lit(2.0) * (lit(1.0) - 
                    (lit(1.0) / (lit(1.0) + exp(lit(-0.717) * _abs(col("wald_statistic")) - lit(0.416) * _pow(_abs(col("wald_statistic")), 2)))))
            )
        )
        
        # Apply Benjamini-Hochberg FDR correction
        df = self.apply_fdr_correction(df, "pvalue", "padj")
        
        return df
    
    def apply_fdr_correction(
        self,
        df: DataFrame,
        pvalue_col: str,
        padj_col: str
    ) -> DataFrame:
        """
        Apply Benjamini-Hochberg FDR correction.
        
        Args:
            df: DataFrame with p-values
            pvalue_col: Name of p-value column
            padj_col: Name for adjusted p-value column
            
        Returns:
            DataFrame with adjusted p-values
        """
        logger.info("Applying Benjamini-Hochberg FDR correction")
        
        # Total number of tests
        n_tests = df.count()
        
        # Rank p-values
        window_spec = Window.orderBy(col(pvalue_col))
        df_ranked = df.withColumn("rank", row_number().over(window_spec))
        
        # Calculate adjusted p-values
        df_ranked = df_ranked.withColumn(
            padj_col,
            when(col(pvalue_col).isNull(), lit(None))
            .otherwise(
                # BH correction: p * n / rank
                when(col(pvalue_col) * lit(n_tests) / col("rank") > 1.0, lit(1.0))
                .otherwise(col(pvalue_col) * lit(n_tests) / col("rank"))
            )
        )
        
        # Enforce monotonicity (cumulative minimum from highest rank)
        window_rev = Window.orderBy(col("rank").desc()).rowsBetween(Window.unboundedPreceding, Window.currentRow)
        df_ranked = df_ranked.withColumn(
            padj_col,
            _min(col(padj_col)).over(window_rev)
        )
        
        return df_ranked.drop("rank")
    
    def filter_significant_genes(
        self,
        df: DataFrame,
        padj_threshold: float = 0.05,
        log2fc_threshold: float = 1.0
    ) -> DataFrame:
        """
        Filter for significantly differentially expressed genes.
        
        Args:
            df: DataFrame with DE results
            padj_threshold: Adjusted p-value threshold
            log2fc_threshold: Absolute log2 fold change threshold
            
        Returns:
            Filtered DataFrame
        """
        logger.info(f"Filtering for padj < {padj_threshold} and |log2FC| > {log2fc_threshold}")
        
        df_sig = df.filter(
            (col("padj") < lit(padj_threshold)) &
            (_abs(col("log2FoldChange")) > lit(log2fc_threshold))
        )
        
        n_sig = df_sig.count()
        logger.info(f"Found {n_sig} significantly differentially expressed genes")
        
        return df_sig
    
    def write_results(self, df: DataFrame, output_path: str) -> None:
        """
        Write results to S3 in Parquet format.
        
        Args:
            df: DataFrame with results
            output_path: S3 path for output
        """
        logger.info(f"Writing results to: {output_path}")
        
        # Select and order columns for output
        output_cols = [
            "gene_id", "gene_name", "chromosome",
            "baseMean", "log2FoldChange", "lfcSE",
            "pvalue", "padj",
            "cond1_mean", "cond2_mean"
        ]
        
        # Only include columns that exist
        existing_cols = [c for c in output_cols if c in df.columns]
        
        df.select(*existing_cols).write.mode("overwrite").parquet(output_path)
        logger.info("Results written successfully")


def main():
    """Main analysis workflow."""
    if len(sys.argv) < 5:
        print("Usage: spark-submit spark_deseq2_analysis.py <input_path> <output_path> <condition1_samples> <condition2_samples>")
        print("Example: spark-submit spark_deseq2_analysis.py s3://bucket/input/ s3://bucket/output/ 'sample_01,sample_02,sample_03' 'sample_04,sample_05,sample_06'")
        sys.exit(1)
    
    input_path = sys.argv[1]
    output_path = sys.argv[2]
    condition1_samples = sys.argv[3].split(',')
    condition2_samples = sys.argv[4].split(',')
    
    logger.info("=== Starting DESeq2-style Differential Expression Analysis ===")
    logger.info(f"Input path: {input_path}")
    logger.info(f"Output path: {output_path}")
    logger.info(f"Condition 1 samples: {condition1_samples}")
    logger.info(f"Condition 2 samples: {condition2_samples}")
    
    # Initialize Spark session
    spark = SparkSession.builder \
        .appName("DESeq2-style Differential Expression Analysis") \
        .config("spark.sql.adaptive.enabled", "true") \
        .config("spark.sql.adaptive.coalescePartitions.enabled", "true") \
        .getOrCreate()
    
    try:
        # Initialize analyzer
        analyzer = DESeq2SparkAnalyzer(spark)
        
        # Read data
        df = analyzer.read_count_data(input_path)
        
        # Get all sample columns
        all_samples = condition1_samples + condition2_samples
        
        # Calculate size factors
        size_factors = analyzer.calculate_size_factors(df, all_samples)
        
        # Normalize counts
        df_norm = analyzer.normalize_counts(df, all_samples, size_factors)
        
        # Calculate log2 fold changes
        df_fc = analyzer.calculate_log2_fold_change(
            df_norm,
            condition1_samples,
            condition2_samples
        )
        
        # Calculate statistics
        df_stats = analyzer.calculate_statistics(
            df_fc,
            condition1_samples,
            condition2_samples
        )
        
        # Add placeholder for lfcSE (would be calculated properly in production)
        df_stats = df_stats.withColumn("lfcSE", lit(0.0))
        
        # Write all results
        analyzer.write_results(df_stats, output_path + "/all_genes")
        
        # Filter and write significant genes
        df_sig = analyzer.filter_significant_genes(df_stats)
        analyzer.write_results(df_sig, output_path + "/significant_genes")
        
        logger.info("=== Analysis completed successfully ===")
        
    except Exception as e:
        logger.error(f"Analysis failed with error: {str(e)}")
        raise
    
    finally:
        spark.stop()


if __name__ == "__main__":
    main()
