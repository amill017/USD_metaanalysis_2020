#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(vegan)
library(mefa)
beta <- as.dist(as.matrix(read.csv(file = args[1], header=TRUE, row.names=1, sep="\t")))
map <- read.csv(file = args[2], header=TRUE, sep="\t")
group <- map$Group
sample <- map$Sample_type
f3 <- paste(group,sample)
pm <- adonis(beta ~ group + sample + group*sample, permutations=999)
pval <- pm[[1]][["Pr(>F)"]]
pd <- betadisper(beta, f3)
rp = vector('expression', 3)
rp[1] = substitute(paste(italic(pvalues)), list(MYOTHERVALUE = format(pval[1], digits = 0)))[2]
rp[2] = substitute(expression(italic(Group) == MYOTHERVALUE), list(MYOTHERVALUE = format(pval[1], digits = 3)))[2]
rp[3] = substitute(expression(italic(Sample_type) == MYOTHERVALUE), list(MYOTHERVALUE = format(pval[2], digits = 3)))[2]
rp[4] = substitute(expression(italic(Group*Sample_type) == MYOTHERVALUE), list(MYOTHERVALUE = format(pval[3], digits = 3)))[2]  
#quartz(width=7)
pdf(file=sprintf("data/downloads/beta_group.pdf"), width=10)
plot(pd, label.cex=1, ellipse=TRUE, hull=FALSE, seg.lty="dashed", cex=0.5, main="Group x Sample Type", cex.main=1.5, cex.sub=0.01, cex.lab=1.1, axes=c(1,2))
legend('topright', legend=rp, bty='n')
legend(x=-0.2, y=0.2, 
  legend = c("USD stool", "Control stool", "USD urine", "Control urine", "USD stone"), 
  col = c("blue", "black", "cyan", "red", "green"), 
  pch = c(4,1,5,2,3), 
  bty = "n", 
  pt.cex = 1.5, 
  cex = 0.8, 
  text.col = "black", 
  horiz = F , 
  inset = c(1, 1))
#quartz.save("data/downloads/beta_group_sampletype.png", type="png")
dev.off()
