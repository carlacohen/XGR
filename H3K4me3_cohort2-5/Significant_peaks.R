#Analysis on H3K4me3-seq cohorts 2-5
#Aim 
#Make supplementary tables of significant DE genes

library(dplyr)
library(stringr)
library(tidyr)

#Read in H3K4me3 data (output from DEseq2)
H3K4me3.CD4 <- read.table("Inputs/H3K4me3.CD4.txt", header = TRUE, row.names = NULL, sep = "\t") 
H3K4me3.CD8 <- read.table("Inputs/H3K4me3.CD8.txt", header = TRUE, row.names = NULL, sep = "\t") 
H3K4me3.CD14 <- read.table("Inputs/H3K4me3.CD14.txt", header = TRUE, row.names = NULL, sep = "\t") 

H3K4me3 <- list(H3K4me3.CD4, H3K4me3.CD8, H3K4me3.CD14)
names(H3K4me3) <- c("CD4", "CD8", "CD14")
head(H3K4me3[1])



#Let's have a look at the data, how many genes fit the conditions pvalue<0.05, FC>1.5?

get_sig_genes <- function (df) {
    df2 <- df %>% 
        filter(log2FoldChange > 0.585 | log2FoldChange < -0.585) %>%
        filter(pvalue < 0.05) 
    return (df2)
}

#Run function on all 3 cell type
H3K4me3.sig <- lapply(H3K4me3, function(x) get_sig_genes(x))
H3K4me3.sig <- lapply (H3K4me3.sig, function (x) rename (x, peak_location = row.names))
head(H3K4me3.sig[[1]])
H3K4me3.CD4.sig <- (H3K4me3.sig[[1]]) 
H3K4me3.CD8.sig <- (H3K4me3.sig[[2]]) 
H3K4me3.CD14.sig <- (H3K4me3.sig[[3]])

#Need to add a column for mean_AS and mean_HV
#Read in limma batch corrected counts
H3K4me3.CD4.counts <- read.table ("./Inputs/limma_batchcorrected_counts/CD4_limma_batchcorrected_counts.txt", 
                              row.names = NULL)
H3K4me3.CD8.counts <- read.table ("./Inputs/limma_batchcorrected_counts/CD8_limma_batchcorrected_counts.txt", 
                              row.names = NULL)
H3K4me3.CD14.counts <- read.table ("./Inputs/limma_batchcorrected_counts/CD14_limma_batchcorrected_counts.txt", 
                               row.names = NULL)
counts <- list(H3K4me3.CD4.counts, H3K4me3.CD8.counts, H3K4me3.CD14.counts)
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

H3K4me3.CD4.counts.means <- as.data.frame (counts_means[1])
H3K4me3.CD8.counts.means <- as.data.frame (counts_means[2])
H3K4me3.CD14.counts.means <- as.data.frame (counts_means[3])
H3K4me3.CD4.sig.means <- H3K4me3.CD4.sig %>% left_join(H3K4me3.CD4.counts.means, by = "peak_location") %>%
    select (peak_location, baseMean, mean_AS, mean_HV, log2FoldChange, lfcSE, stat, pvalue, padj, ProxGene, Dist)

head(H3K4me3.CD4.sig.means)
H3K4me3.CD8.sig.means <- H3K4me3.CD8.sig %>% left_join(H3K4me3.CD8.counts.means, by = "peak_location") %>%
    select (peak_location, baseMean, mean_AS, mean_HV, log2FoldChange, lfcSE, stat, pvalue, padj, ProxGene, Dist)
H3K4me3.CD14.sig.means <- H3K4me3.CD14.sig %>% left_join(H3K4me3.CD14.counts.means, by = "peak_location") %>%
    select (peak_location, baseMean, mean_AS, mean_HV, log2FoldChange, lfcSE, stat, pvalue, padj, ProxGene, Dist)

dir.create("SigDEgenes", showWarnings = FALSE)

write.table (H3K4me3.CD4.sig.means, "./SigDEgenes/H3K4me3_CD4_sigDEgenes.txt", sep = "\t", quote = FALSE,
             col.names = TRUE, row.names = FALSE)
write.table (H3K4me3.CD8.sig.means, "./SigDEgenes/H3K4me3_CD8_sigDEgenes.txt", sep = "\t", quote = FALSE,
             col.names = TRUE, row.names = FALSE)
write.table (H3K4me3.CD14.sig.means, "./SigDEgenes/H3K4me3_CD14_sigDEgenes.txt", sep = "\t", quote = FALSE,
             col.names = TRUE, row.names = FALSE)
