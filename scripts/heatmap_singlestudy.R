#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
options(error=expression(NULL))
library(RColorBrewer)
library(gplots)
data_stool <- read.table(file = args[2], header=TRUE)  
rnames <- data_stool[,1]                            # assign labels in column 1 to "rnames"
mat_data_stool <- data.matrix(data_stool[,2:ncol(data_stool)])  # transform column 2-5 into a matrix
rownames(mat_data_stool) <- rnames ##attach row names to dataset
col_breaks = c(seq(-1,0,length=100),  # for red
   seq(0.001,0.08,length=100),           # for yellow
   seq(0.81,1,length=100))
my_palette <- colorRampPalette(c("red", "yellow", "green", "blue", "violet"))(n = 999)
pdf(file=sprintf("data/downloads/heatmap_stool.pdf"), width=9)
heatmap.2(mat_data_stool, Rowv=FALSE, dendrogram = "none", col=my_palette, trace="none", density.info=c("none"), offsetRow=-45, key.title=NA, lhei=c(1.5,5), margins=c(5,5), adjCol=c(0.5, 1), srtCol=360)
mtext("Taxa", at=0.02, line=-4)
mtext("USD", side=2, at=1.05, line=1.5)
mtext("Healthy", side=4, at=1.05, line=-31)
dev.off()
data_urine <- read.table(file = args[1], header=TRUE)  
rnames <- data_urine[,1]                            # assign labels in column 1 to "rnames"
mat_data_urine <- data.matrix(data_urine[,2:ncol(data_urine)])  # transform column 2-5 into a matrix
rownames(mat_data_urine) <- rnames ##attach row names to dataset
col_breaks = c(seq(-1,0,length=100),  # for red
   seq(0.001,0.08,length=100),           # for yellow
   seq(0.81,1,length=100))
my_palette <- colorRampPalette(c("red", "yellow", "green", "blue", "violet"))(n = 999)
pdf(file=sprintf("data/downloads/heatmap_urine.pdf"), width=9)
heatmap.2(mat_data_urine, Rowv=FALSE, dendrogram = "none", col=my_palette, trace="none", density.info=c("none"), offsetRow=-45, key.title=NA, lhei=c(1.5,5), margins=c(5,5), adjCol=c(0.5, 1), srtCol=360)
mtext("Taxa", at=0.02, line=-4)
mtext("USD", side=2, at=1.05, line=1.5)
mtext("Healthy", side=4, at=1.05, line=-31)
dev.off()

data_stone <- read.table(file = args[3], header=TRUE)  
rnames <- data_stone[,1]                            # assign labels in column 1 to "rnames"
data_stone$Microbiota <- data_stone$Stone
mat_data_stone <- data.matrix(data_stone[,2:ncol(data_stone)])  # transform column 2-5 into a matrix
rownames(mat_data_stone) <- rnames ##attach row names to dataset
col_breaks = c(seq(-1,0,length=100),  # for red
   seq(0.001,0.08,length=100),           # for yellow
   seq(0.81,1,length=100))
my_palette <- colorRampPalette(c("red", "yellow", "green", "blue", "violet"))(n = 999)
pdf(file=sprintf("data/downloads/heatmap_stone.pdf"), width=9)
heatmap.2(mat_data_stone, Rowv=FALSE, dendrogram = "none", col=my_palette, trace="none", density.info=c("none"), offsetRow=-45, key.title=NA, lhei=c(1.5,5), margins=c(5,5), adjCol=c(0.5, 1), srtCol=360)
mtext("Taxa", at=0.02, line=-4)
mtext("Less abundant", side=2, at=1.05, line=1.5)
mtext("More abundant", side=4, at=1.05, line=-31)
dev.off()