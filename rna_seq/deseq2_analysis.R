# DESeq2 Differential Expression Analysis Example
# This script assumes you have run the Python QC/preprocessing and have a counts_for_deseq2.csv file.

# Load required libraries
library(DESeq2)


# Use Snakemake input/output if available
if (exists("snakemake")) {
  in_file <- snakemake@input[[1]]
  out_results <- snakemake@output[[1]]
  out_maplot <- snakemake@output[[2]]
} else {
  in_file <- 'counts_for_deseq2.csv'
  out_results <- 'deseq2_results.csv'
  out_maplot <- 'deseq2_MAplot.png'
}

# Read in count matrix
df <- read.csv(in_file, row.names=1)

# Create a sample table (replace with your actual sample info)
sampleTable <- data.frame(
  row.names = colnames(df),
  condition = c('A', 'A', 'B', 'B')
)

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = df, colData = sampleTable, design = ~ condition)


# Run DESeq2 with fitType='mean' for small datasets
dds <- DESeq(dds, fitType='mean')
res <- results(dds)

# Output results
write.csv(as.data.frame(res), file=out_results)

# Basic MA plot
png(out_maplot)
plotMA(res, main='DESeq2 MA-plot', ylim=c(-2,2))
dev.off()

cat(paste('DESeq2 analysis complete. Results saved to', out_results, '\n'))
