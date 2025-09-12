# PCA Plot in R for bulk RNA-Seq
library(ggplot2)
library(DESeq2)
# Load DESeq2 results or count matrix
dds <- DESeqDataSetFromMatrix(countData = as.matrix(read.csv("../../rna_seq/sample_counts.csv", row.names=1)),
                              colData = DataFrame(condition=rep("A", ncol(read.csv("../../rna_seq/sample_counts.csv"))-1)),
                              design = ~1)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- DESeq(dds)
vsd <- vst(dds)
pcaData <- plotPCA(vsd, intgroup = c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal()
