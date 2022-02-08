##Aim to generate a list of genes that are of future interest
##eg for Felicie's eQTL analysis
setwd("~/R/XGR")
library(dplyr)
#library(stringr)
library(tidyr)
dir.create("Interesting_genes")

###Top DE genes###

#Read in RNA data (output from DEseq2)
RNA.CD4 <- read.table("./RNAseq_cohort2-5/Inputs/RNA.CD4.txt", header = TRUE, row.names = NULL, sep = "\t") 
RNA.CD8 <- read.table("./RNAseq_cohort2-5/Inputs/RNA.CD8.txt", header = TRUE, row.names = NULL, sep = "\t") 
RNA.CD14 <- read.table("./RNAseq_cohort2-5/Inputs/RNA.CD14.txt", header = TRUE, row.names = NULL, sep = "\t") 


#Get list of genes that fit the conditions padj<0.05, FC>1.5?

get_sig_genes <- function (df) {
    df2 <- df %>% 
        filter(log2FoldChange > 0.585 | log2FoldChange < -0.585) %>%
        filter(padj < 0.05) %>%
        select(name)
    return (df2)
}

#Run function on all 3 cell type
RNA.CD4.sig <- get_sig_genes(RNA.CD4)
RNA.CD8.sig <- get_sig_genes(RNA.CD8)
RNA.CD14.sig <- get_sig_genes(RNA.CD14)

write.table(RNA.CD4.sig, "./Interesting_genes/RNA.CD4.sig.txt", quote = FALSE, 
            sep = "\t", col.names = FALSE, row.names = FALSE)
write.table(RNA.CD8.sig, "./Interesting_genes/RNA.CD8.sig.txt", quote = FALSE, 
            sep = "\t", col.names = FALSE, row.names = FALSE)
write.table(RNA.CD14.sig, "./Interesting_genes/RNA.CD14.sig.txt", quote = FALSE, 
            sep = "\t", col.names = FALSE, row.names = FALSE)

###Genes linked to Top200 ATAC peaks###

#Read in ATAC data (output from DEseq2)
ATAC.CD4 <- read.table("./ATAC_cohort2-5/Inputs/ATAC.CD4.txt", header = TRUE, row.names = NULL, sep = "\t") 
ATAC.CD8 <- read.table("./ATAC_cohort2-5/Inputs/ATAC.CD8.txt", header = TRUE, row.names = NULL, sep = "\t") 
ATAC.CD14 <- read.table("./ATAC_cohort2-5/Inputs/ATAC.CD14.txt", header = TRUE, row.names = NULL, sep = "\t") 

#Expand first column to chr start end
ATAC.CD4 <- ATAC.CD4 %>% separate(row.names, c("chr", "start", "end"))
ATAC.CD8 <- ATAC.CD8 %>% separate(row.names, c("chr", "start", "end"))
ATAC.CD14 <- ATAC.CD14 %>% separate(row.names, c("chr", "start", "end"))

#change characters to numeric
ATAC.CD4$start <- as.numeric(ATAC.CD4$start)
ATAC.CD4$end <- as.numeric(ATAC.CD4$end)
ATAC.CD8$start <- as.numeric(ATAC.CD8$start)
ATAC.CD8$end <- as.numeric(ATAC.CD8$end)
ATAC.CD14$start <- as.numeric(ATAC.CD14$start)
ATAC.CD14$end <- as.numeric(ATAC.CD14$end)


Top200peaks <- function (df) {
    df2 <- df %>% 
        filter(log2FoldChange > 0.585 | log2FoldChange < -0.585) %>%
        arrange(padj) %>%
        slice_head(n = 200)
    return (df2)
}

#Run function on all 3 cell type
ATAC.CD4.Top200 <- Top200peaks(ATAC.CD4)
ATAC.CD8.Top200 <- Top200peaks(ATAC.CD8)
ATAC.CD14.Top200 <- Top200peaks(ATAC.CD14)

#Create a list of genes
ATAC.CD4.Top200.proximal <- ATAC.CD4.Top200 %>% select (ProxGene) %>% 
    drop_na() %>% rename("gene" = "ProxGene")
ATAC.CD8.Top200.proximal <- ATAC.CD8.Top200 %>% select (ProxGene) %>% 
    drop_na() %>% rename("gene" = "ProxGene")
ATAC.CD14.Top200.proximal <- ATAC.CD14.Top200 %>% select (ProxGene) %>% 
    drop_na() %>% rename("gene" = "ProxGene")

#Read in data of peaks linked to genes via PCHIC
ATAC.CD4.PCHIC <- read.table("./ATAC_cohort2-5/Inputs/PCHIC/ATAC_ML_CD4_PCHIC.txt", header = FALSE, sep = "\t")
colnames (ATAC.CD4.PCHIC) <- c("chr", "start", "end", "1", "chr_loop", "start_loop", "end_loop", "interaction", "gene", "distance")
ATAC.CD8.PCHIC <- read.table("./ATAC_cohort2-5/Inputs/PCHIC/ATAC_ML_CD8_PCHIC.txt", header = FALSE, sep = "\t")
colnames (ATAC.CD8.PCHIC) <- c("chr", "start", "end", "1", "chr_loop", "start_loop", "end_loop", "interaction", "gene", "distance")
ATAC.CD14.PCHIC <- read.table("./ATAC_cohort2-5/Inputs/PCHIC/ATAC_ML_CD14_PCHIC.txt", header = FALSE, sep = "\t")
colnames (ATAC.CD14.PCHIC) <- c("chr", "start", "end", "1", "chr_loop", "start_loop", "end_loop", "interaction", "gene", "distance")

#Find Top 200 genes that have entries with PCHIC
# Make function to find Top 200 genes that have entries with PCHIC
extract.ATAC.PCHIC <- function(ATAC_Top200peaks, ATAC.ML.PCHIC) {
    ATAC_Top200_PCHIC <- inner_join(ATAC_Top200peaks, ATAC.ML.PCHIC)
    #Create a list of genes for input to XGR
    #Filter out rows where gene is "."
    ATAC_Top200_PCHIC_genes <- ATAC_Top200_PCHIC %>% filter (gene !=".")
    #split genes separated by ";" into a single list
    list <- strsplit(ATAC_Top200_PCHIC_genes$gene, ";")
    # use unlist to put each gene on a new row
    final_list <- as.data.frame(unlist(list))
    #show only unique genes
    unique_list <- unique(final_list)
    #rename columns to "gene"
    colnames(unique_list) <- c("gene")
    return(unique_list)
}
#Run function on 3 cell types
ATAC.CD4.Top200.PCHIC <- extract.ATAC.PCHIC(ATAC.CD4.Top200, ATAC.CD4.PCHIC)
ATAC.CD8.Top200.PCHIC <- extract.ATAC.PCHIC(ATAC.CD8.Top200, ATAC.CD8.PCHIC)
ATAC.CD14.Top200.PCHIC <- extract.ATAC.PCHIC(ATAC.CD14.Top200, ATAC.CD14.PCHIC)
head(ATAC.CD4.Top200.proximal)
head(ATAC.CD4.Top200.PCHIC)

#Make a combined list of ATAC genes linked by proximity and PCHIC
ATAC.CD4.Top200.genes <- rbind (ATAC.CD4.Top200.proximal, ATAC.CD4.Top200.PCHIC)
ATAC.CD8.Top200.genes <- rbind (ATAC.CD8.Top200.proximal, ATAC.CD8.Top200.PCHIC)
ATAC.CD14.Top200.genes <- rbind (ATAC.CD14.Top200.proximal, ATAC.CD14.Top200.PCHIC)

#Save files
write.table(ATAC.CD4.Top200.genes, "./Interesting_genes/ATAC.CD4.Top200.txt", quote = FALSE, 
            sep = "\t", col.names = FALSE, row.names = FALSE)
write.table(ATAC.CD8.Top200.genes, "./Interesting_genes/ATAC.CD8.Top200.txt", quote = FALSE, 
            sep = "\t", col.names = FALSE, row.names = FALSE)
write.table(ATAC.CD14.Top200.genes, "./Interesting_genes/ATAC.CD14.Top200.txt", quote = FALSE, 
            sep = "\t", col.names = FALSE, row.names = FALSE)


###Genes linked to Top200 H3K4me3 peaks###

#Read in H3K4me3 data (output from DEseq2)
H3K4me3.CD4 <- read.table("./H3K4me3_cohort2-5/Inputs/H3K4me3.CD4.txt", header = TRUE, row.names = NULL, sep = "\t") 
H3K4me3.CD8 <- read.table("./H3K4me3_cohort2-5/Inputs/H3K4me3.CD8.txt", header = TRUE, row.names = NULL, sep = "\t") 
H3K4me3.CD14 <- read.table("./H3K4me3_cohort2-5/Inputs/H3K4me3.CD14.txt", header = TRUE, row.names = NULL, sep = "\t") 

#Expand first column to chr start end
H3K4me3.CD4 <- H3K4me3.CD4 %>% separate(row.names, c("chr", "start", "end"))
H3K4me3.CD8 <- H3K4me3.CD8 %>% separate(row.names, c("chr", "start", "end"))
H3K4me3.CD14 <- H3K4me3.CD14 %>% separate(row.names, c("chr", "start", "end"))

#change characters to numeric
H3K4me3.CD4$start <- as.numeric(H3K4me3.CD4$start)
H3K4me3.CD4$end <- as.numeric(H3K4me3.CD4$end)
H3K4me3.CD8$start <- as.numeric(H3K4me3.CD8$start)
H3K4me3.CD8$end <- as.numeric(H3K4me3.CD8$end)
H3K4me3.CD14$start <- as.numeric(H3K4me3.CD14$start)
H3K4me3.CD14$end <- as.numeric(H3K4me3.CD14$end)


Top200peaks <- function (df) {
    df2 <- df %>% 
        filter(log2FoldChange > 0.585 | log2FoldChange < -0.585) %>%
        arrange(padj) %>%
        slice_head(n = 200)
    return (df2)
}

#Run function on all 3 cell type
H3K4me3.CD4.Top200 <- Top200peaks(H3K4me3.CD4)
H3K4me3.CD8.Top200 <- Top200peaks(H3K4me3.CD8)
H3K4me3.CD14.Top200 <- Top200peaks(H3K4me3.CD14)

#Create a list of genes
H3K4me3.CD4.Top200.proximal <- H3K4me3.CD4.Top200 %>% select (ProxGene) %>% 
    drop_na() %>% rename("gene" = "ProxGene")
H3K4me3.CD8.Top200.proximal <- H3K4me3.CD8.Top200 %>% select (ProxGene) %>% 
    drop_na() %>% rename("gene" = "ProxGene")
H3K4me3.CD14.Top200.proximal <- H3K4me3.CD14.Top200 %>% select (ProxGene) %>% 
    drop_na() %>% rename("gene" = "ProxGene")

#Read in data of peaks linked to genes via PCHIC
H3K4me3.CD4.PCHIC <- read.table("./H3K4me3_cohort2-5/Inputs/PCHIC/H3K4me3_ML_CD4_PCHIC.txt", header = FALSE, sep = "\t")
colnames (H3K4me3.CD4.PCHIC) <- c("chr", "start", "end", "1", "chr_loop", "start_loop", "end_loop", "interaction", "gene", "distance")
H3K4me3.CD8.PCHIC <- read.table("./H3K4me3_cohort2-5/Inputs/PCHIC/H3K4me3_ML_CD8_PCHIC.txt", header = FALSE, sep = "\t")
colnames (H3K4me3.CD8.PCHIC) <- c("chr", "start", "end", "1", "chr_loop", "start_loop", "end_loop", "interaction", "gene", "distance")
H3K4me3.CD14.PCHIC <- read.table("./H3K4me3_cohort2-5/Inputs/PCHIC/H3K4me3_ML_CD14_PCHIC.txt", header = FALSE, sep = "\t")
colnames (H3K4me3.CD14.PCHIC) <- c("chr", "start", "end", "1", "chr_loop", "start_loop", "end_loop", "interaction", "gene", "distance")

#Find Top 200 genes that have entries with PCHIC
# Make function to find Top 200 genes that have entries with PCHIC
extract.H3K4me3.PCHIC <- function(H3K4me3_Top200peaks, H3K4me3.ML.PCHIC) {
    H3K4me3_Top200_PCHIC <- inner_join(H3K4me3_Top200peaks, H3K4me3.ML.PCHIC)
    #Create a list of genes for input to XGR
    #Filter out rows where gene is "."
    H3K4me3_Top200_PCHIC_genes <- H3K4me3_Top200_PCHIC %>% filter (gene !=".")
    #split genes separated by ";" into a single list
    list <- strsplit(H3K4me3_Top200_PCHIC_genes$gene, ";")
    # use unlist to put each gene on a new row
    final_list <- as.data.frame(unlist(list))
    #show only unique genes
    unique_list <- unique(final_list)
    #rename columns to "gene"
    colnames(unique_list) <- c("gene")
    return(unique_list)
}
#Run function on 3 cell types
H3K4me3.CD4.Top200.PCHIC <- extract.H3K4me3.PCHIC(H3K4me3.CD4.Top200, H3K4me3.CD4.PCHIC)
H3K4me3.CD8.Top200.PCHIC <- extract.H3K4me3.PCHIC(H3K4me3.CD8.Top200, H3K4me3.CD8.PCHIC)
H3K4me3.CD14.Top200.PCHIC <- extract.H3K4me3.PCHIC(H3K4me3.CD14.Top200, H3K4me3.CD14.PCHIC)
head(H3K4me3.CD4.Top200.proximal)
head(H3K4me3.CD4.Top200.PCHIC)

#Make a combined list of H3K4me3 genes linked by proximity and PCHIC
H3K4me3.CD4.Top200.genes <- rbind (H3K4me3.CD4.Top200.proximal, H3K4me3.CD4.Top200.PCHIC)
H3K4me3.CD8.Top200.genes <- rbind (H3K4me3.CD8.Top200.proximal, H3K4me3.CD8.Top200.PCHIC)
H3K4me3.CD14.Top200.genes <- rbind (H3K4me3.CD14.Top200.proximal, H3K4me3.CD14.Top200.PCHIC)

#Save files
write.table(H3K4me3.CD4.Top200.genes, "./Interesting_genes/H3K4me3.CD4.Top200.txt", quote = FALSE, 
            sep = "\t", col.names = FALSE, row.names = FALSE)
write.table(H3K4me3.CD8.Top200.genes, "./Interesting_genes/H3K4me3.CD8.Top200.txt", quote = FALSE, 
            sep = "\t", col.names = FALSE, row.names = FALSE)
write.table(H3K4me3.CD14.Top200.genes, "./Interesting_genes/H3K4me3.CD14.Top200.txt", quote = FALSE, 
            sep = "\t", col.names = FALSE, row.names = FALSE)

###Genes linked to Top200 H3K27Ac peaks###

#Read in H3K27Ac data (output from DEseq2)
H3K27Ac.CD4 <- read.table("./H3K27Ac_cohort2-5/Inputs/H3K27Ac.CD4.txt", header = TRUE, row.names = NULL, sep = "\t") 
H3K27Ac.CD8 <- read.table("./H3K27Ac_cohort2-5/Inputs/H3K27Ac.CD8.txt", header = TRUE, row.names = NULL, sep = "\t") 
H3K27Ac.CD14 <- read.table("./H3K27Ac_cohort2-5/Inputs/H3K27Ac.CD14.txt", header = TRUE, row.names = NULL, sep = "\t") 

#Expand first column to chr start end
H3K27Ac.CD4 <- H3K27Ac.CD4 %>% separate(row.names, c("chr", "start", "end"))
H3K27Ac.CD8 <- H3K27Ac.CD8 %>% separate(row.names, c("chr", "start", "end"))
H3K27Ac.CD14 <- H3K27Ac.CD14 %>% separate(row.names, c("chr", "start", "end"))

#change characters to numeric
H3K27Ac.CD4$start <- as.numeric(H3K27Ac.CD4$start)
H3K27Ac.CD4$end <- as.numeric(H3K27Ac.CD4$end)
H3K27Ac.CD8$start <- as.numeric(H3K27Ac.CD8$start)
H3K27Ac.CD8$end <- as.numeric(H3K27Ac.CD8$end)
H3K27Ac.CD14$start <- as.numeric(H3K27Ac.CD14$start)
H3K27Ac.CD14$end <- as.numeric(H3K27Ac.CD14$end)


Top200peaks <- function (df) {
    df2 <- df %>% 
        filter(log2FoldChange > 0.585 | log2FoldChange < -0.585) %>%
        arrange(padj) %>%
        slice_head(n = 200)
    return (df2)
}

#Run function on all 3 cell type
H3K27Ac.CD4.Top200 <- Top200peaks(H3K27Ac.CD4)
H3K27Ac.CD8.Top200 <- Top200peaks(H3K27Ac.CD8)
H3K27Ac.CD14.Top200 <- Top200peaks(H3K27Ac.CD14)

#Create a list of genes
H3K27Ac.CD4.Top200.proximal <- H3K27Ac.CD4.Top200 %>% select (ProxGene) %>% 
    drop_na() %>% rename("gene" = "ProxGene")
H3K27Ac.CD8.Top200.proximal <- H3K27Ac.CD8.Top200 %>% select (ProxGene) %>% 
    drop_na() %>% rename("gene" = "ProxGene")
H3K27Ac.CD14.Top200.proximal <- H3K27Ac.CD14.Top200 %>% select (ProxGene) %>% 
    drop_na() %>% rename("gene" = "ProxGene")

#Read in data of peaks linked to genes via PCHIC
H3K27Ac.CD4.PCHIC <- read.table("./H3K27Ac_cohort2-5/Inputs/PCHIC/H3K27Ac_ML_CD4_PCHIC.txt", header = FALSE, sep = "\t")
colnames (H3K27Ac.CD4.PCHIC) <- c("chr", "start", "end", "1", "chr_loop", "start_loop", "end_loop", "interaction", "gene", "distance")
H3K27Ac.CD8.PCHIC <- read.table("./H3K27Ac_cohort2-5/Inputs/PCHIC/H3K27Ac_ML_CD8_PCHIC.txt", header = FALSE, sep = "\t")
colnames (H3K27Ac.CD8.PCHIC) <- c("chr", "start", "end", "1", "chr_loop", "start_loop", "end_loop", "interaction", "gene", "distance")
H3K27Ac.CD14.PCHIC <- read.table("./H3K27Ac_cohort2-5/Inputs/PCHIC/H3K27Ac_ML_CD14_PCHIC.txt", header = FALSE, sep = "\t")
colnames (H3K27Ac.CD14.PCHIC) <- c("chr", "start", "end", "1", "chr_loop", "start_loop", "end_loop", "interaction", "gene", "distance")

#Find Top 200 genes that have entries with PCHIC
# Make function to find Top 200 genes that have entries with PCHIC
extract.H3K27Ac.PCHIC <- function(H3K27Ac_Top200peaks, H3K27Ac.ML.PCHIC) {
    H3K27Ac_Top200_PCHIC <- inner_join(H3K27Ac_Top200peaks, H3K27Ac.ML.PCHIC)
    #Create a list of genes for input to XGR
    #Filter out rows where gene is "."
    H3K27Ac_Top200_PCHIC_genes <- H3K27Ac_Top200_PCHIC %>% filter (gene !=".")
    #split genes separated by ";" into a single list
    list <- strsplit(H3K27Ac_Top200_PCHIC_genes$gene, ";")
    # use unlist to put each gene on a new row
    final_list <- as.data.frame(unlist(list))
    #show only unique genes
    unique_list <- unique(final_list)
    #rename columns to "gene"
    colnames(unique_list) <- c("gene")
    return(unique_list)
}
#Run function on 3 cell types
H3K27Ac.CD4.Top200.PCHIC <- extract.H3K27Ac.PCHIC(H3K27Ac.CD4.Top200, H3K27Ac.CD4.PCHIC)
H3K27Ac.CD8.Top200.PCHIC <- extract.H3K27Ac.PCHIC(H3K27Ac.CD8.Top200, H3K27Ac.CD8.PCHIC)
H3K27Ac.CD14.Top200.PCHIC <- extract.H3K27Ac.PCHIC(H3K27Ac.CD14.Top200, H3K27Ac.CD14.PCHIC)
head(H3K27Ac.CD4.Top200.proximal)
head(H3K27Ac.CD4.Top200.PCHIC)

#Make a combined list of H3K27Ac genes linked by proximity and PCHIC
H3K27Ac.CD4.Top200.genes <- rbind (H3K27Ac.CD4.Top200.proximal, H3K27Ac.CD4.Top200.PCHIC)
H3K27Ac.CD8.Top200.genes <- rbind (H3K27Ac.CD8.Top200.proximal, H3K27Ac.CD8.Top200.PCHIC)
H3K27Ac.CD14.Top200.genes <- rbind (H3K27Ac.CD14.Top200.proximal, H3K27Ac.CD14.Top200.PCHIC)

#Save files
write.table(H3K27Ac.CD4.Top200.genes, "./Interesting_genes/H3K27Ac.CD4.Top200.txt", quote = FALSE, 
            sep = "\t", col.names = FALSE, row.names = FALSE)
write.table(H3K27Ac.CD8.Top200.genes, "./Interesting_genes/H3K27Ac.CD8.Top200.txt", quote = FALSE, 
            sep = "\t", col.names = FALSE, row.names = FALSE)
write.table(H3K27Ac.CD14.Top200.genes, "./Interesting_genes/H3K27Ac.CD14.Top200.txt", quote = FALSE, 
            sep = "\t", col.names = FALSE, row.names = FALSE)

###Genes linked to Top200 eRNA peaks###

#Read in eRNA data (output from DEseq2)
eRNA.CD4 <- read.table("./eRNA_cohort2-5/Inputs/eRNA.CD4.txt", header = TRUE, sep = "\t") 
eRNA.CD8 <- read.table("./eRNA_cohort2-5/Inputs/eRNA.CD8.txt", header = TRUE, sep = "\t") 
eRNA.CD14 <- read.table("./eRNA_cohort2-5/Inputs/eRNA.CD14.txt", header = TRUE, sep = "\t") 

#Expand first column to chr start end
eRNA.CD4 <- eRNA.CD4 %>% separate(peak.name, c("chr", "start", "end"))
eRNA.CD8 <- eRNA.CD8 %>% separate(peak.name, c("chr", "start", "end"))
eRNA.CD14 <- eRNA.CD14 %>% separate(peak.name, c("chr", "start", "end"))

#change characters to numeric
eRNA.CD4$start <- as.numeric(eRNA.CD4$start)
eRNA.CD4$end <- as.numeric(eRNA.CD4$end)
eRNA.CD8$start <- as.numeric(eRNA.CD8$start)
eRNA.CD8$end <- as.numeric(eRNA.CD8$end)
eRNA.CD14$start <- as.numeric(eRNA.CD14$start)
eRNA.CD14$end <- as.numeric(eRNA.CD14$end)



Top200peaks <- function (df) {
    df2 <- df %>% 
        filter(log2FoldChange > 0.585 | log2FoldChange < -0.585) %>%
        arrange(padj) %>%
        slice_head(n = 200)
    return (df2)
}

#Run function on all 3 cell type
eRNA.CD4.Top200 <- Top200peaks(eRNA.CD4)
eRNA.CD8.Top200 <- Top200peaks(eRNA.CD8)
eRNA.CD14.Top200 <- Top200peaks(eRNA.CD14)

#Create a list of genes
eRNA.CD4.Top200.proximal <- eRNA.CD4.Top200 %>% select (ProxGene) %>% 
    drop_na() %>% rename("gene" = "ProxGene")
eRNA.CD8.Top200.proximal <- eRNA.CD8.Top200 %>% select (ProxGene) %>% 
    drop_na() %>% rename("gene" = "ProxGene")
eRNA.CD14.Top200.proximal <- eRNA.CD14.Top200 %>% select (ProxGene) %>% 
    drop_na() %>% rename("gene" = "ProxGene")

#Read in data of peaks linked to genes via PCHIC
eRNA.CD4.PCHIC <- read.table("./eRNA_cohort2-5/Inputs/PCHIC/eRNA_subML_CD4_PCHIC.txt", header = FALSE, sep = "\t")
colnames (eRNA.CD4.PCHIC) <- c("chr", "start", "end", "1", "chr_loop", "start_loop", "end_loop", "interaction", "gene", "distance")
eRNA.CD8.PCHIC <- read.table("./eRNA_cohort2-5/Inputs/PCHIC/eRNA_subML_CD8_PCHIC.txt", header = FALSE, sep = "\t")
colnames (eRNA.CD8.PCHIC) <- c("chr", "start", "end", "1", "chr_loop", "start_loop", "end_loop", "interaction", "gene", "distance")
eRNA.CD14.PCHIC <- read.table("./eRNA_cohort2-5/Inputs/PCHIC/eRNA_subML_CD14_PCHIC.txt", header = FALSE, sep = "\t")
colnames (eRNA.CD14.PCHIC) <- c("chr", "start", "end", "1", "chr_loop", "start_loop", "end_loop", "interaction", "gene", "distance")

#Find Top 200 genes that have entries with PCHIC
# Make function to find Top 200 genes that have entries with PCHIC
extract.eRNA.PCHIC <- function(eRNA_Top200peaks, eRNA.ML.PCHIC) {
    eRNA_Top200_PCHIC <- inner_join(eRNA_Top200peaks, eRNA.ML.PCHIC)
    #Create a list of genes for input to XGR
    #Filter out rows where gene is "."
    eRNA_Top200_PCHIC_genes <- eRNA_Top200_PCHIC %>% filter (gene !=".")
    #split genes separated by ";" into a single list
    list <- strsplit(eRNA_Top200_PCHIC_genes$gene, ";")
    # use unlist to put each gene on a new row
    final_list <- as.data.frame(unlist(list))
    #show only unique genes
    unique_list <- unique(final_list)
    #rename columns to "gene"
    colnames(unique_list) <- c("gene")
    return(unique_list)
}
#Run function on 3 cell types
eRNA.CD4.Top200.PCHIC <- extract.eRNA.PCHIC(eRNA.CD4.Top200, eRNA.CD4.PCHIC)
eRNA.CD8.Top200.PCHIC <- extract.eRNA.PCHIC(eRNA.CD8.Top200, eRNA.CD8.PCHIC)
eRNA.CD14.Top200.PCHIC <- extract.eRNA.PCHIC(eRNA.CD14.Top200, eRNA.CD14.PCHIC)
head(eRNA.CD4.Top200.proximal)
head(eRNA.CD4.Top200.PCHIC)

#Make a combined list of eRNA genes linked by proximity and PCHIC
eRNA.CD4.Top200.genes <- rbind (eRNA.CD4.Top200.proximal, eRNA.CD4.Top200.PCHIC)
eRNA.CD8.Top200.genes <- rbind (eRNA.CD8.Top200.proximal, eRNA.CD8.Top200.PCHIC)
eRNA.CD14.Top200.genes <- rbind (eRNA.CD14.Top200.proximal, eRNA.CD14.Top200.PCHIC)

#Save files
write.table(eRNA.CD4.Top200.genes, "./Interesting_genes/eRNA.CD4.Top200.txt", quote = FALSE, 
            sep = "\t", col.names = FALSE, row.names = FALSE)
write.table(eRNA.CD8.Top200.genes, "./Interesting_genes/eRNA.CD8.Top200.txt", quote = FALSE, 
            sep = "\t", col.names = FALSE, row.names = FALSE)
write.table(eRNA.CD14.Top200.genes, "./Interesting_genes/eRNA.CD14.Top200.txt", quote = FALSE, 
            sep = "\t", col.names = FALSE, row.names = FALSE)

##Alternative splicing###

# read the files
Alt_splicing_sig.CD4 <- read.table("./Alternative_splicing/Cohorts_2-5/Results/Alt_splicing.sig.CD4.txt", header = TRUE, sep = "\t") 
Alt_splicing_sig.CD8 <- read.table("./Alternative_splicing/Cohorts_2-5/Results/Alt_splicing.sig.CD8.txt", header = TRUE, sep = "\t") 
Alt_splicing_sig.CD14 <- read.table("./Alternative_splicing/Cohorts_2-5/Results/Alt_splicing.sig.CD14.txt", header = TRUE, sep = "\t") 

#select the genes
Alt_splicing_sig.CD4.genes <- Alt_splicing_sig.CD4 %>% select (hgnc_symbol) %>%
    drop_na() %>% rename("gene" = "hgnc_symbol") %>% unique()
Alt_splicing_sig.CD8.genes <- Alt_splicing_sig.CD8 %>% select (hgnc_symbol) %>%
    drop_na() %>% rename("gene" = "hgnc_symbol") %>% unique()
Alt_splicing_sig.CD14.genes <- Alt_splicing_sig.CD14 %>% select (hgnc_symbol) %>%
    drop_na() %>% rename("gene" = "hgnc_symbol") %>% unique()

#Save files
write.table(Alt_splicing_sig.CD4.genes, "./Interesting_genes/Alt_splicing.CD4.Top200.txt", quote = FALSE, 
            sep = "\t", col.names = FALSE, row.names = FALSE)
write.table(Alt_splicing_sig.CD8.genes, "./Interesting_genes/Alt_splicing.CD8.Top200.txt", quote = FALSE, 
            sep = "\t", col.names = FALSE, row.names = FALSE)
write.table(Alt_splicing_sig.CD14.genes, "./Interesting_genes/Alt_splicing.CD14.Top200.txt", quote = FALSE, 
            sep = "\t", col.names = FALSE, row.names = FALSE)


###ChromHMM###
#Best approached from individual scripts as fairly complex, see
#./R/XGR/ChromHMM/ChromHMM.CD4/ChromHMM.XGR.PCHIC.proximity.R
#Done only for CD14 as those are the only significant cell type

#Read in results
ChromHMM.promoter.CD14 <- read.table("./Interesting_genes/ChromHMM.promoter.CD14.txt", 
                                     header = FALSE, sep = "\t")

ChromHMM.enhancer.CD14 <- read.table("./Interesting_genes/ChromHMM.enhancer.CD14.txt", 
                                     header = FALSE, sep = "\t")


