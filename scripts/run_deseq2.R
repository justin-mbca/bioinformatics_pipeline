# Differential expression analysis using DESeq2
# Usage: Rscript run_deseq2.R counts.csv sample_info.csv output.csv

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})

args <- commandArgs(trailingOnly=TRUE)
counts_file <- args[1]
sample_info_file <- args[2]
output_file <- args[3]

counts <- read.csv(counts_file, row.names=1)
sample_info <- read.csv(sample_info_file, row.names=1)

dds <- DESeqDataSetFromMatrix(countData=counts, colData=sample_info, design=~condition)
dds <- DESeq(dds)
res <- results(dds)
write.csv(as.data.frame(res), file=output_file)

# PCA plot
pdf("pca.pdf")
vsd <- vst(dds)
plotPCA(vsd, intgroup="condition")
dev.off()

# Volcano plot
pdf("volcano.pdf")
with(as.data.frame(res), plot(log2FoldChange, -log10(padj), pch=20, main="Volcano Plot"))
dev.off()
