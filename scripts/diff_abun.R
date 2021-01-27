#!/usr/bin/env Rscript 
library(DESeq2)
path <- args[5]
args = commandArgs(trailingOnly=TRUE)
cts <- as.matrix(read.csv(args[1], sep="\t", header=TRUE, row.names="OTU_ID"))
storage.mode(cts) = "integer"
map <- read.csv(args[2], sep="\t", row.names=1)
####add 1 to each value in matrix
scts <- cts + 1
dds <- DESeqDataSetFromMatrix(countData= scts, colData=map, design= ~ Group)
dds_res <- DESeq(dds)
results_dds <- results(dds_res)
write.table(results_dds, file=paste(path, "diff_abun_stool.txt", sep = ""))
cts <- as.matrix(read.csv(args[3], sep="\t", header=TRUE, row.names="OTU_ID"))
storage.mode(cts) = "integer"
map <- read.csv(args[4], sep="\t", row.names=1)
####add 1 to each value in matrix
scts <- cts + 1
dds <- DESeqDataSetFromMatrix(countData= scts, colData=map, design= ~ Group)
dds_res <- DESeq(dds)
results_dds <- results(dds_res)
write.table(results_dds, file=paste(path, "diff_abun_urine.txt", sep = ""))