#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
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
x <- as.name(args[5])
r_save <- args[4]
x_save <- args[5]
title <- paste(sample_type,"x", r, "x", x)
samp_save <- args[3]
r3 <- sub_sample[!grepl("blank", sub_sample[[r]]),]
x3 <- r3[!grepl("blank", sub_sample[[x]]),]
r2 <- as.factor(x3[[r]])
x2 <- as.factor(x3[[x]])
r4 <- droplevels(r2)
x4 <- droplevels(x2)
res1 <- aov(PD ~ r4 + x4 + r4*x4, data=x3)
pval <- summary(res1)[[1]][["Pr(>F)"]]
rp = vector('expression', 3)
rp[1] = substitute(paste(italic(pvalues)), list(MYOTHERVALUE = format(pval[1], digits = 0)))[2]
rp[2] = substitute(expression(italic(Var1) == MYOTHERVALUE), list(MYOTHERVALUE = format(pval[1], digits = 3)))[2]
rp[3] = substitute(expression(italic(Var2) == MYOTHERVALUE), list(MYOTHERVALUE = format(pval[2], digits = 3)))[2]
rp[4] = substitute(expression(italic(Var1_Var2) == MYOTHERVALUE), list(MYOTHERVALUE = format(pval[3], digits = 3)))[2]          
#quartz(width=10)
pdf(file=sprintf("data/downloads/alpha_%s_%s_%s.pdf",samp_save,r_save,x_save), width=10)
boxplot(PD~r4*x4, ylab="Phylogenetic Diversity", data=x3, main=c(title), col=(c("red", "blue")))  
leg <- paste(r, rp)
legend('topright', legend=rp, bty='n')
dev.off()
#quartz.save(file=sprintf("data/alpha_%s_%s_%s.png",samp_save,r_save,x_save), type="png")

