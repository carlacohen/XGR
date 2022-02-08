#Use chi squared test to find out if NOTCH pathway genes are statistically more likely
#to be differentially expressed in AS patients

library(tidyverse)

#Read in RNA data (output from DEseq2)
RNA.CD14 <- read.table("Inputs/RNA.CD14.txt", header = TRUE, sep = "\t") 
nrow (RNA.CD14)
#13717 genes are expressed in CD14 monocytes
head (RNA.CD14)
#How many genes are significantly expressed?
RNA.CD14.sig <- RNA.CD14 %>% #filter(log2FoldChange < -0.585 | log2FoldChange > 0.585) %>%
    filter (padj < 0.05) 
nrow(RNA.CD14.sig)
#300 genes are significantly differential

#Import list of Human NOTCH pathway genes

Human_NOTCH_genes <-  read.table("./KEGG_pathway_analysis/HUMAN_NOTCH_genes.txt",header = TRUE, sep = "\t") 
nrow (Human_NOTCH_genes)
#46 genes in NOTCH pathway


#How many NOTCH pathway genes are differentially expressed?
Human_NOTCH_genes
head(RNA.CD14.sig)


RNA.CD14.sig %>% inner_join(Human_NOTCH_genes)
RNA.CD14 %>% inner_join(Human_NOTCH_genes)
