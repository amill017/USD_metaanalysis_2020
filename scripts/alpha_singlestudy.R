#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
alpha <- read.csv(file = args[1], header=TRUE, sep="\t")
map <- read.csv(file = args[2], header=TRUE, sep="\t")
bind <- cbind(alpha, map)
res1 <- aov(PD ~ Group + Sample_type + Group*Sample_type, data=bind)
pval <- summary(res1)[[1]][["Pr(>F)"]]
rp = vector('expression', 3)
rp[1] = substitute(paste(italic(pvalues)), list(MYOTHERVALUE = format(pval[1], digits = 0)))[2]
rp[2] = substitute(expression(italic(Group) == MYOTHERVALUE), list(MYOTHERVALUE = format(pval[1], digits = 3)))[2]
rp[3] = substitute(expression(italic(Sample_type) == MYOTHERVALUE), list(MYOTHERVALUE = format(pval[2], digits = 3)))[2]
rp[4] = substitute(expression(italic(Group*Sample_type) == MYOTHERVALUE), list(MYOTHERVALUE = format(pval[3], digits = 3)))[2]     
pdf(file=sprintf("data/downloads/alpha_group.pdf"),width=10)
#quartz(width=10)
boxplot(PD~Group*Sample_type, ylab="Phylogenetic Diversity", data=bind, col=(c("red", "blue")))  
legend('topright', legend=rp, bty='n')
#quartz.save("data/downloads/alpha.png", type="png")
dev.off()
