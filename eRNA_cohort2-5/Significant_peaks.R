#Analysis on eRNA-seq cohorts 2-5
#Aim 
#Make supplementary tables of significant DE genes

library(dplyr)
library(stringr)
library(tidyr)

#Read in eRNA data (output from DEseq2)
eRNA.CD4 <- read.table("Inputs/eRNA.CD4.txt", header = TRUE, row.names = NULL, sep = "\t") 
eRNA.CD8 <- read.table("Inputs/eRNA.CD8.txt", header = TRUE, row.names = NULL, sep = "\t") 
eRNA.CD14 <- read.table("Inputs/eRNA.CD14.txt", header = TRUE, row.names = NULL, sep = "\t") 

eRNA <- list(eRNA.CD4, eRNA.CD8, eRNA.CD14)
names(eRNA) <- c("CD4", "CD8", "CD14")
head(eRNA[1])



#Let's have a look at the data, how many genes fit the conditions pvalue<0.05, FC>1.5?

get_sig_genes <- function (df) {
    df2 <- df %>% 
        filter(log2FoldChange > 0.585 | log2FoldChange < -0.585) %>%
        filter(pvalue < 0.05) 
    return (df2)
}

#Run function on all 3 cell type
eRNA.sig <- lapply(eRNA, function(x) get_sig_genes(x))
#eRNA.sig <- lapply (eRNA.sig, function (x) rename (x, peak_location = row.names))
head(eRNA.sig[[1]])
eRNA.CD4.sig <- (eRNA.sig[[1]]) 
eRNA.CD8.sig <- (eRNA.sig[[2]]) 
eRNA.CD14.sig <- (eRNA.sig[[3]])

#Need to add a column for mean_AS and mean_HV
#Read in limma batch corrected counts
eRNA.CD4.counts <- read.table ("./Inputs/limma_batchcorrected_counts/CD4_limma_batchcorrected_counts.txt", 
                              row.names = NULL)
eRNA.CD8.counts <- read.table ("./Inputs/limma_batchcorrected_counts/CD8_limma_batchcorrected_counts.txt", 
                              row.names = NULL)
eRNA.CD14.counts <- read.table ("./Inputs/limma_batchcorrected_counts/CD14_limma_batchcorrected_counts.txt", 
                               row.names = NULL)
counts <- list(eRNA.CD4.counts, eRNA.CD8.counts, eRNA.CD14.counts)
#names(counts) <- c("CD4", "CD8", "CD14")
head(counts[1])


get_means <- function (df) {
    df2 <- df %>% 
    mutate(mean_AS = rowMeans(select(., starts_with("AS")))) %>%
    mutate(mean_HV = rowMeans(select(., starts_with("X"), starts_with("HV"), starts_with("CTL"))))
}

counts_means <- lapply(counts, function(x) get_means(x))
counts_means <- lapply (counts_means, function (x) rename (x, peak.name = row.names))
head(counts_means[1])

eRNA.CD4.counts.means <- as.data.frame (counts_means[1])
eRNA.CD8.counts.means <- as.data.frame (counts_means[2])
eRNA.CD14.counts.means <- as.data.frame (counts_means[3])
eRNA.CD4.sig.means <- eRNA.CD4.sig %>% left_join(eRNA.CD4.counts.means, by = "peak.name") %>%
    select (peak.name, baseMean, mean_AS, mean_HV, log2FoldChange, lfcSE, stat, pvalue, padj, ProxGene, Dist)

head(eRNA.CD4.sig.means)
eRNA.CD8.sig.means <- eRNA.CD8.sig %>% left_join(eRNA.CD8.counts.means, by = "peak.name") %>%
    select (peak.name, baseMean, mean_AS, mean_HV, log2FoldChange, lfcSE, stat, pvalue, padj, ProxGene, Dist)
eRNA.CD14.sig.means <- eRNA.CD14.sig %>% left_join(eRNA.CD14.counts.means, by = "peak.name") %>%
    select (peak.name, baseMean, mean_AS, mean_HV, log2FoldChange, lfcSE, stat, pvalue, padj, ProxGene, Dist)

dir.create("SigDEgenes", showWarnings = FALSE)

write.table (eRNA.CD4.sig.means, "./SigDEgenes/eRNA_CD4_sigDEgenes.txt", sep = "\t", quote = FALSE,
             col.names = TRUE, row.names = FALSE)
write.table (eRNA.CD8.sig.means, "./SigDEgenes/eRNA_CD8_sigDEgenes.txt", sep = "\t", quote = FALSE,
             col.names = TRUE, row.names = FALSE)
write.table (eRNA.CD14.sig.means, "./SigDEgenes/eRNA_CD14_sigDEgenes.txt", sep = "\t", quote = FALSE,
             col.names = TRUE, row.names = FALSE)
