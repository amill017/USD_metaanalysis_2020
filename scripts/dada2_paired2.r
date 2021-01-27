#!/usr/bin/env Rscript
##dada2 pipeline
args = commandArgs(trailingOnly=TRUE)
options(error=expression(NULL))
library(dada2)
##Put path to directory of reads in variable
path <- args[1]
##Separate forward and reverse reads - change this for already paired reads
fnFs <- sort(list.files(path, pattern="_R1.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
names(filtFs) <- sample.names
out <- filterAndTrim(fnFs, filtFs, truncLen=c(240),
              maxN=0, maxEE=c(10), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
##Learn error rates              
errF <- learnErrors(filtFs, multithread=TRUE)
derepFs <- derepFastq(filtFs, verbose=TRUE)
names(derepFs) <- sample.names
##Sample inference
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
##Make count table
seqtab <- makeSequenceTable(dadaFs)
##Remove Chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
##Assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "new_refseqs.fna", multithread=TRUE)
##Hand off to phyloseq
library(phyloseq)
library(Biostrings)
library(ggplot2)
library(microbiome)
library(metagMisc)
library(plyr)
library(DECIPHER)
samples.out <- rownames(seqtab.nochim)
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),  
               tax_table(taxa))
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
##Filter chlorophasts and mitochondria
ps <- subset_taxa(ps, (Class!="Chloroplast") | is.na(Class))
ps <- subset_taxa(ps, (Order!="Rickettsiales") | is.na(Order))
ps <- subset_taxa(ps, (Kingdom!="Eukaryota") | is.na(Kingdom))
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
##Generate phylogenetic tree
seqs <- getSequences(seqtab)
names(seqs) <- seqs
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)
library(msa)
library(phangorn)
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)

## negative edges length changed to 0!

fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                       rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)
pst <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),  
               tax_table(taxa), phy_tree(fitGTR$tree))
pst <- subset_taxa(pst, (Class!="Chloroplast") | is.na(Class))
pst <- subset_taxa(pst, (Order!="Rickettsiales") | is.na(Order))
pst <- subset_taxa(pst, (Kingdom!="Eukaryota") | is.na(Kingdom))
taxa_names(pst) <- paste0("ASV", seq(ntaxa(ps)))
tree1 = phy_tree(fitGTR$tree)
taxa_names(tree1) <- paste0("ASV", seq(ntaxa(ps)))
ape::write.tree(tree1, file=paste(path, "/rep_set.tre", sep = ""))
##Write asv tables to file
write_phyloseq(pst, type = "OTU", path = args[1])
write_phyloseq(pst, type = "TAXONOMY", path = args[1])
file.rename(paste(path, "/otu_table.csv", sep = ""), paste(path, "/otu_table_raw.csv", sep = ""))
file.rename(paste(path, "/taxonomy_table.csv", sep = ""), paste(path, "/taxonomy_table_raw.csv", sep = ""))
##log2 transform data
ps_log2 <- transform_sample_counts(pst, log2)
otu_table(ps_log2)[otu_table(ps_log2) < 0.0] <- 0.0
##Write asv tables to file
write_phyloseq(ps_log2, type = "OTU", path = args[1])
write_phyloseq(ps_log2, type = "TAXONOMY", path = args[1])
file.rename(paste(path, "/otu_table.csv", sep = ""), paste(path, "/otu_table_normal.csv", sep = ""))
file.rename(paste(path, "/taxonomy_table.csv", sep = ""), paste(path, "/taxonomy_table_normal.csv", sep = ""))
library(PhyloMeasures)
pd_alpha <- phyloseq_phylo_div(ps_log2, measures = c("PD"))
write.table(pd_alpha, file=paste(path, "/current_alpha.txt", sep = ""))
##Unifrac analysis
wunifrac <- as.data.frame(as.matrix(phyloseq::distance(ps_log2, method="wunifrac")))
write.table(wunifrac, file=paste(path, "/wunifrac.txt", sep = ""))