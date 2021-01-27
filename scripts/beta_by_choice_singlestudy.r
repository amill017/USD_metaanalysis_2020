#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
##Beta diversity by choice
##Run as follows
##Rscript beta_by_choice.r [sample_type (stool, urine, stone, or all)] [metadata group 1] [metadata group 2]
library(vegan)
library(mefa)
beta <- as.dist(read.csv(args[1], header=TRUE, row.names=1, sep="\t"))
map <- read.csv(args[2], header=TRUE, sep="\t", stringsAsFactors = FALSE)
args <- commandArgs(TRUE)
sample_type <- as.name(args[3])
r <- as.name(args[4])
x <- as.name(args[5])
type_save <- args[3]
r_save <- args[4]
x_save <- args[5]
title <- paste(sample_type,"x", r, "x", x)
if(sample_type == "stool"){
sub_sample <- map[grepl("stool", map$Sample_type),]
} else if(sample_type == "urine"){
sub_sample <- map[grepl("urine", map$Sample_type),]
} else if(sample_type == "stone"){
sub_sample <- map[grepl("stone", map$Sample_type),]
} else {
sub_sample <- map
}
r2 <- sub_sample[[r]]
x2 <- sub_sample[[x]]
Sample <- sub_sample$X.SampleID
bind <- as.data.frame(cbind(Sample, r2, x2))
r4 <- bind[!grepl("blank", bind$r2),]  ##dm filter
r5 <- r4[!grepl("blank", r4$x2),]  ##dm filter 
##Subset beta
subs <- as.factor(r5$Sample)  ##if choose urine, all are urine
subs2 <- droplevels(subs) ##if choose urine, all are urine
subs3 <- as.vector(subs2)
library(usedist)
beta_sub <- dist_subset(beta, subs3)  ##Currently, if choose urine, all are feces
r2 <- as.factor(r5$r2)
x2 <- as.factor(r5$x2)
r4 <- droplevels(r2)
x4 <- droplevels(x2)
varbind <- data.frame(r4,x4)
f3 <- paste(r4,x4)
pm_var1 <- adonis(beta_sub ~ r4, permutations=999)
pm_var2 <- adonis(beta_sub ~ x4, permutations=999)
pm_combined <- adonis(beta_sub ~ r4*x4, permutations=999) 
pval1 <- pm_var1[[1]][["Pr(>F)"]]
pval2 <- pm_var2[[1]][["Pr(>F)"]]
pval3 <- pm_combined[[1]][["Pr(>F)"]]
pval32 <- pval3[[3]]
pd <- betadisper(beta_sub, f3)
rp = vector('expression', 3)
rp[1] = substitute(expression(italic(Var1) == MYOTHERVALUE), list(MYOTHERVALUE = format(pval1, digits = 3)))[2]
rp[2] = substitute(expression(italic(Var2) == MYOTHERVALUE), list(MYOTHERVALUE = format(pval2, digits = 3)))[2]
rp[3] = substitute(expression(italic(Var1_Var2) == MYOTHERVALUE), list(MYOTHERVALUE = format(pval32, digits = 3)))[2]
#quartz()
pdf(file=sprintf("data/downloads/beta_%s_%s_%s.pdf",type_save,r_save,x_save))
plot(pd, label.cex=0.95, ellipse=TRUE, hull=FALSE, seg.lty="dashed", cex=0.5, main=c(title), cex.main=1.5, cex.sub=0.01, cex.lab=1.1, axes=c(1,2))
mtext("p-values", at=0.9)
legend('topright', legend=rp, bty='n')
#quartz.save(file=sprintf("full_data/beta_%s_%s_%s.png",type_save,r_save,x_save), type="png")
dev.off()

