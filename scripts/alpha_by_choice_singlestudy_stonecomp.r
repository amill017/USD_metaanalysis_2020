#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(ggplot2)
##Alpha diversity by choice
##Run as follows
##Rscript alpha_by_choice.r [sample_type (stool, urine, stone, or all)] [metadata group 1] [metadata group 2]
alpha <- read.csv(args[1], header=TRUE, sep="\t")
map <- read.csv(args[2], header=TRUE, sep="\t")
bind <- cbind(alpha, map)
args <- commandArgs(TRUE)
##Select sample type
sample_type <- as.name(args[3])
if(sample_type == "stool"){
sub_sample <- bind[grepl("stool", bind$Sample_type),]
} else if(sample_type == "urine"){
sub_sample <- bind[grepl("urine", bind$Sample_type),]
} else if(sample_type == "stone"){
sub_sample <- bind[grepl("stone", bind$Sample_type),]
} else {
sub_sample <- bind
}
r <- as.name(args[4])
#x <- as.name(args[5])
r_save <- args[4]
#x_save <- args[5]
#title <- paste(sample_type,"x", r)
samp_save <- args[3]
r3 <- sub_sample[!grepl("blank", sub_sample[[r]]),]
#x3 <- r3[!grepl("blank", sub_sample[[x]]),]
r2 <- as.factor(r3[[r]])
#x2 <- as.factor(x3[[x]])
r4 <- droplevels(r2)
#x4 <- droplevels(x2)
res1 <- aov(PD ~ r4, data=r3)
pval <- summary(res1)[[1]][["Pr(>F)"]]
rp = vector('expression', 1)
rp[1] = substitute(paste(italic(pvalues)), list(MYOTHERVALUE = format(pval[1], digits = 0)))[2]
rp[2] = substitute(expression(italic(Var1) == MYOTHERVALUE), list(MYOTHERVALUE = format(pval[1], digits = 3)))[2]
#rp[3] = substitute(expression(italic(Var2) == MYOTHERVALUE), list(MYOTHERVALUE = format(pval[2], digits = 3)))[2]
#rp[4] = substitute(expression(italic(Var1_Var2) == MYOTHERVALUE), list(MYOTHERVALUE = format(pval[3], digits = 3)))[2]
#quartz(width=10)
pdf(file=sprintf("data/downloads/alpha_%s_%s_wpval.pdf",samp_save,r_save), width=6)
boxplot(PD~r4, ylab="Phylogenetic Diversity", data=r3, col=(c("red", "blue")))
leg <- paste(r, rp)
legend('topright', legend=rp, bty='n')
dev.off()
pdf(file=sprintf("data/downloads/alpha_%s_%s.pdf",samp_save,r_save), width=6)
ggplot(r3, aes(y = PD, x = r4, fill = r4, color = r4)) +
  geom_boxplot(alpha = 1, width=0.75, color="black", fill="white", outlier.shape=NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1)) + xlab(r) + labs(color=r) +ylab("Phylogenetic Diversity") +
    theme(text = element_text(size=20))
  #dev.off()
#dev.off()
#quartz.save(file=sprintf("data/alpha_%s_%s_%s.png",samp_save,r_save), type="png")
