#Analysis on ATAC-seq cohorts 2-5
#Aim 
#Make supplementary tables of significant DE genes

library(dplyr)
library(stringr)
library(tidyr)

#Read in ATAC data (output from DEseq2)
ATAC.CD4 <- read.table("Inputs/ATAC.CD4.txt", header = TRUE, row.names = NULL, sep = "\t") 
ATAC.CD8 <- read.table("Inputs/ATAC.CD8.txt", header = TRUE, row.names = NULL, sep = "\t") 
ATAC.CD14 <- read.table("Inputs/ATAC.CD14.txt", header = TRUE, row.names = NULL, sep = "\t") 

ATAC <- list(ATAC.CD4, ATAC.CD8, ATAC.CD14)
names(ATAC) <- c("CD4", "CD8", "CD14")
head(ATAC[1])



#Let's have a look at the data, how many genes fit the conditions pvalue<0.05, FC>1.5?

get_sig_genes <- function (df) {
    df2 <- df %>% 
        filter(log2FoldChange > 0.585 | log2FoldChange < -0.585) %>%
        filter(pvalue < 0.05) 
    return (df2)
}

#Run function on all 3 cell type
ATAC.sig <- lapply(ATAC, function(x) get_sig_genes(x))
ATAC.sig <- lapply (ATAC.sig, function (x) rename (x, peak_location = row.names))
head(ATAC.sig[[1]])
ATAC.CD4.sig <- (ATAC.sig[[1]]) 
ATAC.CD8.sig <- (ATAC.sig[[2]]) 
ATAC.CD14.sig <- (ATAC.sig[[3]])

#Need to add a column for mean_AS and mean_HV
#Read in limma batch corrected counts
ATAC.CD4.counts <- read.table ("./Inputs/limma_batchcorrected_counts/CD4_limma_batchcorrected_counts.txt", 
                              row.names = NULL)
ATAC.CD8.counts <- read.table ("./Inputs/limma_batchcorrected_counts/CD8_limma_batchcorrected_counts.txt", 
                              row.names = NULL)
ATAC.CD14.counts <- read.table ("./Inputs/limma_batchcorrected_counts/CD14_limma_batchcorrected_counts.txt", 
                               row.names = NULL)
counts <- list(ATAC.CD4.counts, ATAC.CD8.counts, ATAC.CD14.counts)
#names(counts) <- c("CD4", "CD8", "CD14")
head(counts[1])


get_means <- function (df) {
    df2 <- df %>% 
    mutate(mean_AS = rowMeans(select(., starts_with("AS")))) %>%
    mutate(mean_HV = rowMeans(select(., starts_with("X"), starts_with("HV"), starts_with("CTL"))))
}

counts_means <- lapply(counts, function(x) get_means(x))
counts_means <- lapply (counts_means, function (x) rename (x, peak_location = row.names))
head(counts_means[1])

ATAC.CD4.counts.means <- as.data.frame (counts_means[1])
ATAC.CD8.counts.means <- as.data.frame (counts_means[2])
ATAC.CD14.counts.means <- as.data.frame (counts_means[3])
ATAC.CD4.sig.means <- ATAC.CD4.sig %>% left_join(ATAC.CD4.counts.means, by = "peak_location") %>%
    select (peak_location, baseMean, mean_AS, mean_HV, log2FoldChange, lfcSE, stat, pvalue, padj, ProxGene, Dist)

head(ATAC.CD4.sig.means)
ATAC.CD8.sig.means <- ATAC.CD8.sig %>% left_join(ATAC.CD8.counts.means, by = "peak_location") %>%
    select (peak_location, baseMean, mean_AS, mean_HV, log2FoldChange, lfcSE, stat, pvalue, padj, ProxGene, Dist)
ATAC.CD14.sig.means <- ATAC.CD14.sig %>% left_join(ATAC.CD14.counts.means, by = "peak_location") %>%
    select (peak_location, baseMean, mean_AS, mean_HV, log2FoldChange, lfcSE, stat, pvalue, padj, ProxGene, Dist)

dir.create("SigDEgenes", showWarnings = FALSE)

write.table (ATAC.CD4.sig.means, "./SigDEgenes/ATAC_CD4_sigDEgenes.txt", sep = "\t", quote = FALSE,
             col.names = TRUE, row.names = FALSE)
write.table (ATAC.CD8.sig.means, "./SigDEgenes/ATAC_CD8_sigDEgenes.txt", sep = "\t", quote = FALSE,
             col.names = TRUE, row.names = FALSE)
write.table (ATAC.CD14.sig.means, "./SigDEgenes/ATAC_CD14_sigDEgenes.txt", sep = "\t", quote = FALSE,
             col.names = TRUE, row.names = FALSE)
