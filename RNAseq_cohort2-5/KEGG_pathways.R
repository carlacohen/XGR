#Script to make input for KEGG colour mapper to get KEGG pathways coloured by log2Fold change

library(tidyverse)
library(dplyr)

#Read in RNA data (output from DEseq2)
RNA.CD4 <- read.table("Inputs/RNA.CD4.txt", header = TRUE, sep = "\t") 
RNA.CD8 <- read.table("Inputs/RNA.CD8.txt", header = TRUE, sep = "\t") 
RNA.CD14 <- read.table("Inputs/RNA.CD14.txt", header = TRUE, sep = "\t") 

#Read in genes from TNF-a pathway from https://www.genome.jp/dbget-bin/www_bget?pathway:hsa04668
#TNF_genes<- read.table("./KEGG_pathway_analysis/TNF-a_pathway_genes.txt", 
#                         header = FALSE, sep = "\t")
#Read in genes from NOTCH pathway from https://www.genome.jp/dbget-bin/www_bget?pathway:hsa04330
#NOTCH_genes <- read.table("./KEGG_pathway_analysis/NOTCH_pathway_genes.txt", 
#                          header = FALSE, sep = "\t")
#Read in genes from CYTOKINE-CYTOKINE RECEPTOR INTERACTION pathway  hsa04060
#from https://www.kegg.jp/dbget-bin/www_bget?hsa04060
CYTOKINE_genes <- read.table("./KEGG_pathway_analysis/CYTOKINE_pathway_genes.txt", 
                             header = FALSE, sep = "\t")

#Create some fake genes to set the scale consistently
name <- c("Gene1", "Gene2")
CYTOKINE_log2FoldChange <- c(-2, 4)
CYTOKINE_Fake_genes <- data.frame(name, CYTOKINE_log2FoldChange)
colnames(CYTOKINE_Fake_genes) <-  c("name", "log2FoldChange")

pathway_function <- function (pathway, cell_type){
    pathway_genes <- get(paste0(pathway, "_genes"))
    RNA.df <- get(paste0("RNA.", cell_type))
    Fake_genes <- get(paste0(pathway, "_Fake_genes"))
    print(head (pathway_genes))
    print(Fake_genes)
    #format columns of pathway genes df
    colnames(pathway_genes) <- c("KEGG_id", "name")
    pathway_genes.df <- pathway_genes %>% 
        mutate(name =sub(pattern = "(^.*)\\;.*", replacement = "\\1", pathway_genes$name ))%>% 
        select(KEGG_id, name)
    print(head(pathway_genes.df))
    #Find genes in common between RNA DEseq2 output and TNF pathway
    join.df <- pathway_genes.df %>% inner_join(RNA.df) %>%
        select(KEGG_id, log2FoldChange) %>%
        mutate(KEGG_id = paste0("hsa:",KEGG_id)) %>%
        arrange(log2FoldChange)
    print(head(join.df))
    #save the output
    write.table(join.df, paste0("./KEGG_pathway_analysis/", pathway, "_", cell_type, ".txt"), 
                row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE) 
    #make another df for plotting
    names <- pathway_genes.df %>% inner_join(RNA.df) %>%
        select(name, log2FoldChange) %>% rbind(Fake_genes)
    #plot to get correct colour for each gene
    plot <- ggplot (names, aes (x = name, y = log2FoldChange, colour = log2FoldChange)) +
        geom_point(size = 5)+
        scale_colour_gradient2(low = "#00008b", high = "#8b0000")+
        coord_flip()
    ggsave(paste0(pathway, "_", cell_type, "_plot.pdf"), 
           width = 4, height = 16, useDingbats = FALSE, path = "KEGG_pathway_analysis")
    return(plot)
}

pathway_function("CYTOKINE", "CD4")
pathway_function("CYTOKINE", "CD8")
pathway_function("CYTOKINE", "CD14")

#Use this output on the KEGG web browser tool
#https://www.genome.jp/kegg/tool/map_pathway3.html
#Can't work out how to get file upload to work correctly
#Copy and paste gene list scores with a header eg. #hsa TNF

###Need to incorporate different gene list for human pathway
##Have created a new file "HUMAN_NOTCH_genes.txt" generated from human pathway
#Requires an amended function since input is a list of gene names, not number and gene description
#NB this includes some alternate gene names, see Human_NOTCH_genes_alternative_names.xls
Human_NOTCH_genes <-  read.table("./KEGG_pathway_analysis/HUMAN_NOTCH_genes.txt",header = TRUE, sep = "\t") 


Human_NOTCH_log2FoldChange <- c(-1, 1)
Human_NOTCH_Fake_genes <- data.frame(name, Human_NOTCH_log2FoldChange)
colnames(Human_NOTCH_Fake_genes) <-  c("name", "log2FoldChange")

#Sublist for CXC pathway 
#NB this includes some alternate gene names, see CXC_genes_alternative_names.xls
CXC_genes <- read.table("./KEGG_pathway_analysis/CXCR_subfamily_genes.txt", header = TRUE, sep = "\t") 

CXC_log2FoldChange <- c(-2, 4)
CXC_Fake_genes <- data.frame(name, CXC_log2FoldChange)
colnames(CXC_Fake_genes) <-  c("name", "log2FoldChange")

pathway_function_2 <- function (pathway, cell_type){
    pathway_genes <- get(paste0(pathway, "_genes"))
    RNA.df <- get(paste0("RNA.", cell_type))
    Fake_genes <- get(paste0(pathway, "_Fake_genes"))
    print(head (pathway_genes))
    print(Fake_genes)
    #Find genes in common between RNA DEseq2 output and TNF pathway
    join.df <- pathway_genes %>% inner_join(RNA.df) %>%
        select(name, log2FoldChange) %>%
        arrange(log2FoldChange)
    print(head(join.df))
    #save the output
    write.table(join.df, paste0("./KEGG_pathway_analysis/", pathway, "_", cell_type, ".txt"), 
                row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE) 
    #make another df for plotting
    names <- pathway_genes %>% inner_join(RNA.df) %>%
        select(name, log2FoldChange) %>% rbind(Fake_genes)
    #plot to get correct colour for each gene
    plot <- ggplot (names, aes (x = name, y = log2FoldChange, colour = log2FoldChange)) +
        geom_point(size = 5)+
        scale_colour_gradient2(low = "#00008b", high = "#8b0000")+
        coord_flip()
    ggsave(paste0(pathway, "_", cell_type, "_plot.pdf"), 
           width = 4, height = 8, useDingbats = FALSE, path = "KEGG_pathway_analysis")
    return(plot)
}

pathway_function_2("Human_NOTCH", "CD4")
pathway_function_2("Human_NOTCH", "CD8")
pathway_function_2("Human_NOTCH", "CD14")
pathway_function_2("CXC", "CD4")
pathway_function_2("CXC", "CD8")
pathway_function_2("CXC", "CD14")
