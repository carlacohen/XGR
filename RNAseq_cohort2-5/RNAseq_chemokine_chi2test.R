#Use chi squared test to find out if NOTCH pathway genes are statistically more likely
#to be differentially expressed in AS patients

library(tidyverse)

#Read in RNA data (output from DEseq2)
RNA.CD4 <- read.table("Inputs/RNA.CD4.txt", header = TRUE, sep = "\t") 
RNA.CD8 <- read.table("Inputs/RNA.CD8.txt", header = TRUE, sep = "\t") 
RNA.CD14 <- read.table("Inputs/RNA.CD14.txt", header = TRUE, sep = "\t") 

nrow (RNA.CD4) #14496 genes are expressed in CD4
nrow (RNA.CD8) #14423 genes are expressed in CD8
nrow (RNA.CD14)#13717 genes are expressed in CD14 monocytes

#How many genes are significantly expressed?
RNA.CD4.sig <- RNA.CD4 %>% filter(log2FoldChange < -0.585 | log2FoldChange > 0.585) %>%
    filter (padj < 0.05) 
nrow(RNA.CD4.sig) #122 genes are significantly DE in CD14RNA.CD14.sig <- RNA.CD14 %>% filter(log2FoldChange < -0.585 | log2FoldChange > 0.585) %>%
RNA.CD8.sig <- RNA.CD8 %>% filter(log2FoldChange < -0.585 | log2FoldChange > 0.585) %>%
    filter (padj < 0.05) 
nrow(RNA.CD8.sig) #299 genes are significantly DE in CD14RNA.CD14.sig <- RNA.CD14 %>% filter(log2FoldChange < -0.585 | log2FoldChange > 0.585) %>%
RNA.CD14.sig <- RNA.CD14 %>% filter(log2FoldChange < -0.585 | log2FoldChange > 0.585) %>%
    filter (padj < 0.05) 
nrow(RNA.CD14.sig) #300 genes are significantly DE in CD14

#Import list of Human CYTOKINE pathway genes

CYTOKINE_genes <- read.table("./KEGG_pathway_analysis/CYTOKINE_pathway_genes.txt", 
                             header = FALSE, sep = "\t") 
#Add headers and alter formatting of second column
colnames(CYTOKINE_genes) <- c("KEGG_id", "name")
CYTOKINE_genes <- CYTOKINE_genes %>% 
    mutate(name =sub(pattern = "(^.*)\\;.*", replacement = "\\1", CYTOKINE_genes$name ))%>% 
    select(name)

nrow (CYTOKINE_genes)#295 genes in CYTOKINE pathway


#How many CYTOKINE pathway genes are differentially expressed?
RNA.CD4.sig %>% inner_join(CYTOKINE_genes) #9 genes are DE in CD14
RNA.CD8.sig %>% inner_join(CYTOKINE_genes) #14 genes are DE in CD14
RNA.CD14.sig %>% inner_join(CYTOKINE_genes) #10 genes are DE in CD14

#Import list of CXC family genes

CXC_genes <- read.table("./KEGG_pathway_analysis/CXCR_subfamily_genes.txt", header = TRUE, sep = "\t") 
#25 CXC genes

RNA.CD4.sig %>% inner_join(CXC_genes) #3 genes are DE in CD14
RNA.CD8.sig %>% inner_join(CXC_genes) #5 genes are DE in CD14
RNA.CD14.sig %>% inner_join(CXC_genes) #1 genes are DE in CD14


