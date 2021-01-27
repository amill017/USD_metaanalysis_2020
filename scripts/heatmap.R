  ##1/11/18
library(RColorBrewer)
library(gplots)
data <- read.table("data/scaled_difference_urine5_sorted2.txt", header=TRUE)  
rnames <- data[,1]                            # assign labels in column 1 to "rnames"
mat_data <- data.matrix(data[,2:ncol(data)])  # transform column 2-5 into a matrix
rownames(mat_data) <- rnames ##attach row names to dataset
col_breaks = c(seq(-1,0,length=100),  # for red
   seq(0.001,0.08,length=100),           # for yellow
   seq(0.81,1,length=100))
my_palette <- colorRampPalette(c("red", "yellow", "green", "blue", "violet"))(n = 999)
pdf(file=sprintf("data/heatmap_urine.pdf"), width=9)
heatmap.2(mat_data, Rowv=FALSE, dendrogram = "none", col=my_palette, trace="none", density.info=c("none"), offsetRow=-45, key.title=NA, lhei=c(1.5,5), margins=c(5,5), adjCol=c(0.5, 1), srtCol=360)
mtext("Taxa", at=0.02, line=-4)
mtext("USD", side=2, at=1.05, line=1.5)
mtext("Healthy", side=4, at=1.05, line=-31)
dev.off()
data <- read.table("data/scaled_difference_stool5_sorted2.txt", header=TRUE)  
rnames <- data[,1]                            # assign labels in column 1 to "rnames"
mat_data <- data.matrix(data[,2:ncol(data)])  # transform column 2-5 into a matrix
rownames(mat_data) <- rnames ##attach row names to dataset
col_breaks = c(seq(-1,0,length=100),  # for red
   seq(0.001,0.08,length=100),           # for yellow
   seq(0.81,1,length=100))
my_palette <- colorRampPalette(c("red", "yellow", "green", "blue", "violet"))(n = 999)
pdf(file=sprintf("data/heatmap_stool.pdf"), width=9)
heatmap.2(mat_data, Rowv=FALSE, dendrogram = "none", col=my_palette, trace="none", density.info=c("none"), offsetRow=-45, key.title=NA, lhei=c(1.5,5), margins=c(5,5), adjCol=c(0.5, 1), srtCol=360)
mtext("Taxa", at=0.02, line=-4)
mtext("USD", side=2, at=1.05, line=1.5)
mtext("Healthy", side=4, at=1.05, line=-31)
dev.off()
data <- read.table("data/stone_relabun_heatmap3.txt", header=TRUE)
rnames <- data[,1]                            # assign labels in column 1 to "rnames"
mat_data <- data.matrix(data[,2:ncol(data)])  # transform column 2-5 into a matrix
rownames(mat_data) <- rnames ##attach row names to dataset
col_breaks = c(seq(-1,0,length=100),  # for red
   seq(0.001,0.08,length=100),           # for yellow
   seq(0.81,1,length=100))
my_palette <- colorRampPalette(c("red", "yellow", "green", "blue", "violet"))(n = 999)
mat_data <- cbind(mat_data, mat_data)
pdf(file=sprintf("data/heatmap_stone.pdf"), width=9)
heatmap.2(mat_data, Rowv=FALSE, dendrogram = "none", col=my_palette, trace="none", density.info=c("none"), offsetRow=-45, key.title=NA, lhei=c(1.5,5), margins=c(5,5), adjCol=c(0.5, 1), srtCol=360, labCol="")
mtext("OTU", at=0.02, line=-4)
mtext("Less Abundant", side=2, at=1.01, line=1.5)
mtext("More Abundant", side=4, at=1.01, line=-31)
dev.off()
