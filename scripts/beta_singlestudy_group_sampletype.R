#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(vegan)
library(mefa)
beta <- as.dist(read.csv(file = args[1], header=TRUE, row.names=1, sep="\t"))
map <- read.csv(file = args[2], header=TRUE, sep="\t")
#group <- map$Group
sample <- map$Sample_type
#F3 <- paste(group, sample)
pm <- adonis(beta ~ sample, permutations=999)
pval <- pm[[1]][["Pr(>F)"]]
pd <- betadisper(beta, sample)
#pm_var1 <- adonis(beta ~ group, permutations=999)
pm_var2 <- adonis(beta ~ sample, permutations=999)
#pm_combined <- adonis(beta ~ group*sample, permutations=999)
#pval1 <- pm_var1[[1]][["Pr(>F)"]]
pval2 <- pm_var2[[1]][["Pr(>F)"]]
#pval3 <- pm_combined[[1]][["Pr(>F)"]]
#pval32 <- pval3[[3]]
rp = vector('expression', 3)
#rp[1] = substitute(expression(italic(Group) == MYOTHERVALUE), list(MYOTHERVALUE = format(pval1, digits = 3)))[2]
rp[2] = substitute(expression(italic(Sample_type) == MYOTHERVALUE), list(MYOTHERVALUE = format(pval2, digits = 3)))[2]
#rp[3] = substitute(expression(italic(Group*Sample_type) == MYOTHERVALUE), list(MYOTHERVALUE = format(pval32, digits = 3)))[2]
#quartz(width=7)
pdf(file=sprintf("data/downloads/beta_group_sampletype.pdf"), width=10)
plot(pd, label.cex=1.5, ellipse=TRUE, hull=FALSE, seg.lty="dashed", cex=0.5, main="Group x Sample type", cex.main=1.5, cex.sub=0.01, cex.lab=1.1, axes=c(1,2))
mtext("p-values", at=0.5)
legend('topright', legend=rp, bty='n')
legend(x=-0.9, y=1,
  legend = c("USD", "Control"),
  col = c("blue", "black"),
  pch = c(4,1,5,2,3),
  bty = "n",
  pt.cex = 1.5,
  cex = 0.8,
  text.col = "black",
  horiz = F ,
  inset = c(0.1, 0.1))
#quartz.save("data/downloads/beta_group.png", type="png")
dev.off()
