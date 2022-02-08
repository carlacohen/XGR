#Analysis on RNA-seq cohorts 2-5
#Aim 
#Find which DE genes are near SNPs

setwd("~/R/XGR/RNAseq_cohort2-5")
library(dplyr)

#Read in RNA data (output from DEseq2)
RNA.CD4 <- read.table("Inputs/RNA.CD4.txt", header = TRUE, row.names = NULL, sep = "\t") 
RNA.CD8 <- read.table("Inputs/RNA.CD8.txt", header = TRUE, row.names = NULL, sep = "\t") 
RNA.CD14 <- read.table("Inputs/RNA.CD14.txt", header = TRUE, row.names = NULL, sep = "\t") 

RNA <- list(RNA.CD4, RNA.CD8, RNA.CD14)
names(RNA) <- c("CD4", "CD8", "CD14")
head(RNA[1])

#How many genes fit the conditions padj<0.05, FC>1.5?

get_sig_genes <- function (df) {
    df2 <- df %>% 
        filter(log2FoldChange > 0.585 | log2FoldChange < -0.585) %>%
        filter(padj < 0.05) 
    return (df2)
}

#Run function on all 3 cell type
RNA.sig <- lapply(RNA, function(x) get_sig_genes(x))
RNA.sig <- lapply (RNA.sig, function (x) rename (x, EnsemblID = row.names))
head(RNA.sig[[1]])
RNA.CD4.sig <- (RNA.sig[[1]]) #122 genes
RNA.CD8.sig <- (RNA.sig[[2]]) #299 genes
RNA.CD14.sig <- (RNA.sig[[3]]) #300 genes

#Read in locations of genes within 500kb of lead SNP (found using bedtools)
#See U:/AS_analysis/RNAseq_cohorts2-5/DEgenes_near_SNPs

Genes.SNPs.500kb <- read.table ("~/R/Gene_locations/hg19_genes_vs_SNPs.500000.txt")

colnames(Genes.SNPs.500kb) <- c("SNP.chr", "SNP.start", "SNP.end", "SNP", "Gene.chr", "Gene.start", "Gene.end", "name")

#Which genes are within 500kb of SNP and also differentially expressed?

Genes.SNPs.500kb.CD4 <- Genes.SNPs.500kb %>% inner_join(RNA.CD4.sig) %>%
    #add column called "distance" from gene to SNP
    mutate(distance = Gene.start - SNP.start)
write.table(Genes.SNPs.500kb.CD4, "./gene_locations/RNA.CD4.sig.near_SNP.txt", sep = "\t", quote = FALSE, row.names = FALSE)
Genes.SNPs.500kb.CD8 <- Genes.SNPs.500kb %>% inner_join(RNA.CD8.sig) %>%
    mutate(distance = Gene.start - SNP.start)
write.table(Genes.SNPs.500kb.CD8, "./gene_locations/RNA.CD8.sig.near_SNP.txt", sep = "\t", quote = FALSE, row.names = FALSE)
Genes.SNPs.500kb.CD14 <- Genes.SNPs.500kb %>% inner_join(RNA.CD14.sig) %>%
    mutate(distance = Gene.start - SNP.start) 
write.table(Genes.SNPs.500kb.CD14, "./gene_locations/RNA.CD14.sig.near_SNP.txt", sep = "\t", quote = FALSE, row.names = FALSE)


#Overall how many genes are there regardless of cell type?
Allgenes.SNPs.500kb <- rbind(Genes.SNPs.500kb.CD4, Genes.SNPs.500kb.CD8, Genes.SNPs.500kb.CD14) 
length(unique(Allgenes.SNPs.500kb$name))
#35 unique DE genes