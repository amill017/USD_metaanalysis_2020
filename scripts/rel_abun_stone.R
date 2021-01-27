#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
options(error=expression(NULL))
stone <- read.table(args[1], header=TRUE, row.names=1)
stone$taxonomy <- NULL
results <- sapply(stone, function(i) i/sum(i))
write.table(results, args[2], sep="\t", row.names=F)