#Analysis on RNA-seq cohorts 2-5
#Aim 
#Make supplementary tables of significant DE genes


library(dplyr)
library(stringr)
library(tidyr)
#library(GenomicRanges)
library(XGR)

#Define colours
Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000") 


#Read in RNA data (output from DEseq2)
RNA.CD4 <- read.table("Inputs/RNA.CD4.txt", header = TRUE, row.names = NULL, sep = "\t") 
RNA.CD8 <- read.table("Inputs/RNA.CD8.txt", header = TRUE, row.names = NULL, sep = "\t") 
RNA.CD14 <- read.table("Inputs/RNA.CD14.txt", header = TRUE, row.names = NULL, sep = "\t") 

RNA <- list(RNA.CD4, RNA.CD8, RNA.CD14)
names(RNA) <- c("CD4", "CD8", "CD14")
head(RNA[1])

#Let's have a look at the data, how many genes fit the conditions padj<0.05, FC>1.5?

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

#Need to add a column for mean_AS and mean_HV
#Read in limma batch corrected counts
RNA.CD4.counts <- read.table ("./Inputs/limma_batchcorrected_counts/CD4_limma_batchcorrected_counts.txt", 
                              row.names = NULL)
RNA.CD8.counts <- read.table ("./Inputs/limma_batchcorrected_counts/CD8_limma_batchcorrected_counts.txt", 
                              row.names = NULL)
RNA.CD14.counts <- read.table ("./Inputs/limma_batchcorrected_counts/CD14_limma_batchcorrected_counts.txt", 
                               row.names = NULL)
counts <- list(RNA.CD4.counts, RNA.CD8.counts, RNA.CD14.counts)
#names(counts) <- c("CD4", "CD8", "CD14")
head(counts[1])


get_means <- function (df) {
    df2 <- df %>% 
    mutate(mean_AS = rowMeans(select(., starts_with("AS")))) %>%
    mutate(mean_HV = rowMeans(select(., starts_with("X"), starts_with("HV"), starts_with("CTL"))))
}

counts_means <- lapply(counts, function(x) get_means(x))
counts_means <- lapply (counts_means, function (x) rename (x, EnsemblID = row.names))
head(counts_means[1])

RNA.CD4.counts.means <- as.data.frame (counts_means[1])
RNA.CD8.counts.means <- as.data.frame (counts_means[2])
RNA.CD14.counts.means <- as.data.frame (counts_means[3])
RNA.CD4.sig.means <- RNA.CD4.sig %>% left_join(RNA.CD4.counts.means, by = "EnsemblID") %>%
    select (EnsemblID, baseMean, mean_AS, mean_HV, log2FoldChange, lfcSE, stat, pvalue, padj, name)
RNA.CD8.sig.means <- RNA.CD8.sig %>% left_join(RNA.CD8.counts.means, by = "EnsemblID") %>%
    select (EnsemblID, baseMean, mean_AS, mean_HV, log2FoldChange, lfcSE, stat, pvalue, padj, name)
RNA.CD14.sig.means <- RNA.CD14.sig %>% left_join(RNA.CD14.counts.means, by = "EnsemblID") %>%
    select (EnsemblID, baseMean, mean_AS, mean_HV, log2FoldChange, lfcSE, stat, pvalue, padj, name)

dir.create("SigDEgenes", showWarnings = FALSE)

write.table (RNA.CD4.sig.means, "./SigDEgenes/RNA_CD4_sigDEgenes.txt", sep = "\t", quote = FALSE,
             col.names = TRUE, row.names = FALSE)
write.table (RNA.CD8.sig.means, "./SigDEgenes/RNA_CD8_sigDEgenes.txt", sep = "\t", quote = FALSE,
             col.names = TRUE, row.names = FALSE)
write.table (RNA.CD14.sig.means, "./SigDEgenes/RNA_CD14_sigDEgenes.txt", sep = "\t", quote = FALSE,
             col.names = TRUE, row.names = FALSE)
