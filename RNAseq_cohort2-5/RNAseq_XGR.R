#Analysis on RNA-seq cohorts 2-5

setwd("~/R/XGR/RNAseq_cohort2-5")
#Aim 
# Run XGR on DE genes

library(dplyr)
library(stringr)
library(tidyr)
#library(GenomicRanges)
library(XGR)

#Define colours
Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000") 


#Read in RNA data (output from DEseq2)
RNA.CD4 <- read.table("Inputs/RNA.CD4.txt", header = TRUE, sep = "\t") 
RNA.CD8 <- read.table("Inputs/RNA.CD8.txt", header = TRUE, sep = "\t") 
RNA.CD14 <- read.table("Inputs/RNA.CD14.txt", header = TRUE, sep = "\t") 

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
head(RNA.sig[[1]])
nrow(RNA.sig[[1]]) #122 genes
nrow(RNA.sig[[2]]) #299 genes
nrow(RNA.sig[[3]]) #300 genes
RNA.CD4.sig <- as.data.frame(RNA.sig[1]) 
colnames(RNA.CD4.sig)
colnames(RNA.CD4) <- colnames(RNA.CD4.sig)
RNA.CD8.sig <- as.data.frame(RNA.sig[2])
colnames(RNA.CD8.sig)
colnames(RNA.CD8) <- colnames(RNA.CD8.sig)
RNA.CD14.sig <- as.data.frame(RNA.sig[3])
colnames(RNA.CD14) <- colnames(RNA.CD14.sig)
#Therefore perform GO analyis on significant genes, not Top200 genes. 

#Find genes that are significant and 
#(i) common to all 3 cell types
#(ii) unique to each cell type

common.sig.genes <- RNA.CD4.sig %>% inner_join(RNA.CD8.sig, by = c("CD4.name" = "CD8.name")) %>%
    inner_join(RNA.CD14.sig, by = c("CD4.name" = "CD14.name")) %>%
    select (CD4.name, CD4.log2FoldChange, CD4.padj, CD8.log2FoldChange, CD8.padj, CD14.log2FoldChange, CD14.padj)

unique_genes_CD4 <- RNA.CD4.sig %>% full_join(RNA.CD8, by = c("CD4.name" = "CD8.name")) %>%
    full_join(RNA.CD14, by = c("CD4.name" = "CD14.name")) %>%
    select (CD4.name, CD4.baseMean, CD4.log2FoldChange, CD4.padj, 
           CD8.padj, 
          CD14.padj) %>%
    filter (CD4.padj < 0.05) %>%
    filter (CD4.log2FoldChange > 0.585 | CD4.log2FoldChange < -0.585) %>%
    filter (CD8.padj > 0.05) %>%
    filter (CD14.padj > 0.05)
#e.g. TMEM30B is only upregulated in CD4 


#Get list of Top200 genes to compare with cohorts 1-3
#Create a function to select top 200 genes with FC > 1.5 and sorted by padj value
Top200peaks <- function (df) {
    df2 <- df %>% 
        filter(log2FoldChange > 0.585 | log2FoldChange < -0.585) %>%
        arrange(padj) %>%
        slice_head(n = 200)
    return (df2)
}

#Run function on all 3 cell type
RNA.Top200 <- lapply(RNA, function(x) Top200peaks(x))

#Get gene lists
RNA.Top200.genes <- lapply (RNA.Top200, function (x) dplyr::select(x, name))
names(RNA.Top200.genes)

dir.create("Top200_genes")
for (i in 1:length(RNA.Top200.genes)){
    write.table (RNA.Top200.genes[[i]], 
                 paste0("./Top200_genes/",names(RNA.sig.genes.bed[i]),".txt"),
                 col.names = FALSE, 
                 row.names = FALSE,
                 quote = FALSE,
                 sep = "\t")
}


#Background lists of genes 
RNA.background <- lapply (RNA, function (x) dplyr::select (x, name))
RNA.background 
names(RNA.background)

#Convert background gene lists to vectors
RNA.CD4.background <- unlist(RNA.background[[1]])
RNA.CD8.background <- unlist(RNA.background[[2]])
RNA.CD14.background <- unlist(RNA.background[[3]])

#Make a concatenated background set using rbind
RNA.background.concat <- rbind (RNA.background [[1]], 
                            RNA.background [[2]],
                            RNA.background [[3]])
class(RNA.background)
head(RNA.background)

#Convert gene lists to vectors
RNA.CD4.sig <- unlist(RNA.sig[[1]])
RNA.CD8.sig <- unlist(RNA.sig[[2]])
RNA.CD14.sig <- unlist(RNA.sig[[3]])

#Generate background for all 3 cell types
RNA.background <- as.vector(RNA.background.concat$name)
class(RNA.background)

library(XGR)
dir.create("GOBP")
dir.create("Reactome")
dir.create("KEGG")
RData.location <- "http://galahad.well.ox.ac.uk/bigdata"

XGR_function <- function(genelist, background, ontology) {
    #Run xEnricher function 
    #NB settings: for <50 genes use size.range = c(5,2000) min.overlap = 3
    #For larger gene sets than 50 use size.range = c(10,2000) min.overlap = 5 (which is the default setting)
    enriched_genes <- xEnricherGenes(data = genelist, 
                                     background = background, 
                                     test = "hypergeo", 
                                     ontology = ontology, 
                                     RData.location = RData.location,
                                     size.range = c(10,2000),
                                     min.overlap = 5)
    enriched_genes_concise <- xEnrichConciser(enriched_genes)
    return (enriched_genes_concise)
}


RNA.CD4.sig.GOBP <- XGR_function(RNA.CD4.sig, RNA.background, ontology = "GOBP")
RNA.CD4.sig.GOBP.table <- xEnrichViewer(RNA.CD4.sig.GOBP, details = TRUE, 
                                           top_num = 50)

RNA.CD4.sig.GOBP.barplot <- xEnrichBarplot(RNA.CD4.sig.GOBP)
RNA.CD4.sig.GOBP.barplot
write.table(RNA.CD4.sig.GOBP.table, "./GOBP/RNA.CD4.sig.GOBP.txt", sep = "\t", quote = FALSE)

RNA.CD8.sig.GOBP <- XGR_function(RNA.CD8.sig, RNA.background, ontology = "GOBP")
RNA.CD8.sig.GOBP.table <- xEnrichViewer(RNA.CD8.sig.GOBP, details = TRUE, 
                                           top_num = 50)
RNA.CD8.sig.GOBP.table
write.table(RNA.CD8.sig.GOBP.table, "./GOBP/RNA.CD8.sig.GOBP.txt", sep = "\t", quote = FALSE)
RNA.CD8.sig.GOBP.barplot <- xEnrichBarplot(RNA.CD8.sig.GOBP, displayBy="fdr", bar.color = "lightblue-blue",)
RNA.CD8.sig.GOBP.barplot

RNA.CD14.sig.GOBP <- XGR_function(RNA.CD14.sig, RNA.background, ontology = "GOBP")
RNA.CD14.sig.GOBP.table <- xEnrichViewer(RNA.CD14.sig.GOBP, details = TRUE, 
                                     top_num = 50)
write.table(RNA.CD14.sig.GOBP.table, "./GOBP/RNA.CD14.sig.GOBP.txt", sep = "\t", quote = FALSE)

#Compare 3 cell types
RNA.sig.GOBP <- list(RNA.CD4.sig.GOBP, RNA.CD8.sig.GOBP, RNA.CD14.sig.GOBP)
names(RNA.sig.GOBP) <- c("CD4+ T cells", "CD8+ T cells", "CD14+ monocytes")
#Set the order
RNA.sig.GOBP <- RNA.sig.GOBP[c("CD14+ monocytes", "CD8+ T cells", "CD4+ T cells")]
RNA.GOBP.sig.compare <- xEnrichCompare(RNA.sig.GOBP, FDR.cutoff = 0.01, 
                                           signature = FALSE, 
                                           bar.label = FALSE, 
                                           displayBy = "adjp")
RNA.GOBP.sig.compare <- 
    RNA.GOBP.sig.compare + theme(axis.text.y=element_text(size=12,color="black")) +
    scale_fill_manual(values = c(Okabe_Ito[3], Okabe_Ito[5], Okabe_Ito[6])) 
RNA.GOBP.sig.compare

dir.create("figures", showWarnings = FALSE)
ggsave("RNA.sig.GOBP.pdf", RNA.GOBP.sig.compare, width = 12, height = 14, useDingbats = FALSE, path = "./figures/")
ggsave("RNA.sig.GOBP.png", RNA.GOBP.sig.compare, width = 12, height = 14, device = "png", path = "./figures/")


#Make heatmap of results 
#Read in results tables so don't have to re-run all above code
RNA.CD4.sig.GOBP.table <- read.table ("./GOBP/RNA.CD4.sig.GOBP.txt", 
                                   header = TRUE, sep = "\t")
RNA.CD8.sig.GOBP.table <- read.table ("./GOBP/RNA.CD8.sig.GOBP.txt", 
                                   header = TRUE, sep = "\t")
RNA.CD14.GOBP.sig.table <- read.table ("./GOBP/RNA.CD14.sig.GOBP.txt", 
                                    header = TRUE, sep = "\t")

#First collect all the bits of data into a df then convert to a matrix
RNA.CD4.sig.GOBP.padj <- data.frame(RNA.CD4.sig.GOBP.table$name, 
                                 RNA.CD4.sig.GOBP.table$adjp)
colnames(RNA.CD4.sig.GOBP.padj) <- c("name", "CD4.padj")
RNA.CD4.sig.GOBP.padj <- RNA.CD4.sig.GOBP.padj %>% filter(CD4.padj < 0.01)

RNA.CD8.sig.GOBP.padj <- data.frame(RNA.CD8.sig.GOBP.table$name, 
                                       RNA.CD8.sig.GOBP.table$adjp)
colnames(RNA.CD8.sig.GOBP.padj) <- c("name", "CD8.padj")
RNA.CD8.sig.GOBP.padj <- RNA.CD8.sig.GOBP.padj %>% filter(CD8.padj < 0.01)

RNA.CD14.sig.GOBP.padj <- data.frame(RNA.CD14.sig.GOBP.table$name, 
                                       RNA.CD14.sig.GOBP.table$adjp)
colnames(RNA.CD14.sig.GOBP.padj) <- c("name", "CD14.padj")
RNA.CD14.sig.GOBP.padj <- RNA.CD14.sig.GOBP.padj %>% filter(CD14.padj < 0.01)

RNA.sig.GOBP.padj <- RNA.CD4.sig.GOBP.padj %>%
    full_join (RNA.CD8.sig.GOBP.padj) %>%
    full_join (RNA.CD14.sig.GOBP.padj) %>%
    #Add columns for -log10 padj
    mutate(CD4 = -log10(CD4.padj)) %>%
    mutate(CD8 = -log10(CD8.padj)) %>%
    mutate(CD14 = -log10(CD14.padj)) %>%
    select (c(1, 5, 6, 7))

head(RNA.sig.GOBP.padj)
#Table of -log10(padj) --> make heatmap in GraphPad Prism!
write.table(RNA.sig.GOBP.padj, 
            "./GOBP/RNA.sig.GOBP.padj.txt", 
            quote = FALSE, 
            sep = "\t", 
            row.names = FALSE, 
            na = "")



#Now do for Reactome ontology = "MsigdbC2REACTOME"

RNA.CD4.sig.Reactome <- XGR_function(RNA.CD4.sig, RNA.background, ontology = "MsigdbC2REACTOME")
RNA.CD4.sig.Reactome.table <- xEnrichViewer(RNA.CD4.sig.Reactome, details = TRUE, 
                                        top_num = 50)

RNA.CD4.sig.Reactome.barplot <- xEnrichBarplot(RNA.CD4.sig.Reactome)
RNA.CD4.sig.Reactome.barplot

write.table(RNA.CD4.sig.Reactome.table, "./Reactome/RNA.CD4.sig.Reactome.txt", sep = "\t", quote = FALSE)

RNA.CD8.sig.Reactome <- XGR_function(RNA.CD8.sig, RNA.background, ontology = "MsigdbC2REACTOME")
RNA.CD8.sig.Reactome.table <- xEnrichViewer(RNA.CD8.sig.Reactome, details = TRUE, 
                                        top_num = 50)
RNA.CD8.sig.Reactome.table
write.table(RNA.CD8.sig.Reactome.table, "./Reactome/RNA.CD8.sig.Reactome.txt", sep = "\t", quote = FALSE)
RNA.CD8.sig.Reactome.barplot <- xEnrichBarplot(RNA.CD8.sig.Reactome, displayBy="fdr", bar.color = "lightblue-blue",)
RNA.CD8.sig.Reactome.barplot

RNA.CD14.sig.Reactome <- XGR_function(RNA.CD14.sig, RNA.background, ontology = "MsigdbC2REACTOME")
RNA.CD14.sig.Reactome.table <- xEnrichViewer(RNA.CD14.sig.Reactome, details = TRUE, 
                                         top_num = 50)
write.table(RNA.CD14.sig.Reactome.table, "./Reactome/RNA.CD14.sig.Reactome.txt", sep = "\t", quote = FALSE)

RNA.sig.Reactome <- list(RNA.CD4.sig.Reactome, RNA.CD8.sig.Reactome, RNA.CD14.sig.Reactome)
names(RNA.sig.Reactome) <- c("CD4+ T cells", "CD8+ T cells", "CD14+ monocytes")
#Set the order
RNA.sig.Reactome <- RNA.sig.Reactome[c("CD14+ monocytes", "CD8+ T cells", "CD4+ T cells")]
RNA.Reactome.sig.compare <- xEnrichCompare(RNA.sig.Reactome, FDR.cutoff = 0.01, 
                                           signature = FALSE, 
                                           bar.label = FALSE, 
                                           displayBy = "adjp")
RNA.Reactome.sig.compare <- 
    RNA.Reactome.sig.compare + theme(axis.text.y=element_text(size=18,color="black")) +
    scale_fill_manual(values = c(Okabe_Ito[3], Okabe_Ito[5], Okabe_Ito[6])) 
RNA.Reactome.sig.compare

dir.create("figures", showWarnings = FALSE)
ggsave("RNA.sig.Reactome.pdf", RNA.Reactome.sig.compare, width = 12, height = 10, useDingbats = FALSE, path = "./figures/")
ggsave("RNA.sig.Reactome.png", RNA.Reactome.sig.compare, width = 12, height = 10, device = "png", path = "./figures/")



#Make heatmap of results for Top200
#Read in results tables so don't have to re-run all above code

#Make heatmap of results for sig
#Read in results tables so don't have to re-run all above code
RNA.CD4.sig.Reactome.table <- read.table ("./Reactome/RNA.CD4.sig.Reactome.txt", 
                                      header = TRUE, sep = "\t")
RNA.CD8.sig.Reactome.table <- read.table ("./Reactome/RNA.CD8.sig.Reactome.txt", 
                                      header = TRUE, sep = "\t")
RNA.CD14.Reactome.sig.table <- read.table ("./Reactome/RNA.CD14.sig.Reactome.txt", 
                                       header = TRUE, sep = "\t")

#First collect all the bits of data into a df then convert to a matrix
RNA.CD4.sig.Reactome.padj <- data.frame(RNA.CD4.sig.Reactome.table$name, 
                                    RNA.CD4.sig.Reactome.table$adjp)
colnames(RNA.CD4.sig.Reactome.padj) <- c("name", "CD4.padj")
RNA.CD4.sig.Reactome.padj <- RNA.CD4.sig.Reactome.padj %>% filter(CD4.padj < 0.05)

RNA.CD8.sig.Reactome.padj <- data.frame(RNA.CD8.sig.Reactome.table$name, 
                                    RNA.CD8.sig.Reactome.table$adjp)
colnames(RNA.CD8.sig.Reactome.padj) <- c("name", "CD8.padj")
RNA.CD8.sig.Reactome.padj <- RNA.CD8.sig.Reactome.padj %>% filter(CD8.padj < 0.05)

RNA.CD14.sig.Reactome.padj <- data.frame(RNA.CD14.sig.Reactome.table$name, 
                                     RNA.CD14.sig.Reactome.table$adjp)
colnames(RNA.CD14.sig.Reactome.padj) <- c("name", "CD14.padj")
RNA.CD14.sig.Reactome.padj <- RNA.CD14.sig.Reactome.padj %>% filter(CD14.padj < 0.05)

RNA.sig.Reactome.padj <- RNA.CD4.sig.Reactome.padj %>%
    full_join (RNA.CD8.sig.Reactome.padj) %>%
    full_join (RNA.CD14.sig.Reactome.padj) %>%
    #Add columns for -log10 padj
    mutate(CD4 = -log10(CD4.padj)) %>%
    mutate(CD8 = -log10(CD8.padj)) %>%
    mutate(CD14 = -log10(CD14.padj)) %>%
    select (c(1, 5, 6, 7))

RNA.sig.Reactome.padj
#Table of -log10(padj) --> make heatmap in GraphPad Prism!
write.table(RNA.sig.Reactome.padj, 
            "./Reactome/RNA.sig.Reactome.padj.txt", 
            quote = FALSE, 
            sep = "\t", 
            row.names = FALSE, 
            na = "")


#Do for KEGG pathways - ontology = "MsigdbC2KEGG"
RNA.CD4.sig.KEGG <- XGR_function(RNA.CD4.sig, RNA.background, ontology = "MsigdbC2KEGG")
RNA.CD4.sig.KEGG.table <- xEnrichViewer(RNA.CD4.sig.KEGG, details = TRUE, 
                                        top_num = 50)
write.table(RNA.CD4.sig.KEGG.table, "./KEGG/RNA.CD4.sig.KEGG.txt", sep = "\t", quote = FALSE)

RNA.CD8.sig.KEGG <- XGR_function(RNA.CD8.sig, RNA.background, ontology = "MsigdbC2KEGG")
RNA.CD8.sig.KEGG.table <- xEnrichViewer(RNA.CD8.sig.KEGG, details = TRUE, 
                                        top_num = 50)
write.table(RNA.CD8.sig.KEGG.table, "./KEGG/RNA.CD8.sig.KEGG.txt", sep = "\t", quote = FALSE)

RNA.CD14.sig.KEGG <- XGR_function(RNA.CD14.sig, RNA.background, ontology = "MsigdbC2KEGG")
RNA.CD14.sig.KEGG.table <- xEnrichViewer(RNA.CD14.sig.KEGG, details = TRUE, 
                                         top_num = 50)
write.table(RNA.CD14.sig.KEGG.table, "./KEGG/RNA.CD14.sig.KEGG.txt", sep = "\t", quote = FALSE)

#Compare 3 cell types
RNA.sig.KEGG <- list(RNA.CD4.sig.KEGG, RNA.CD8.sig.KEGG, RNA.CD14.sig.KEGG)
names(RNA.sig.KEGG) <- c("CD4+ T cells", "CD8+ T cells", "CD14+ monocytes")
#Set the order
RNA.sig.KEGG <- RNA.sig.KEGG[c("CD14+ monocytes", "CD8+ T cells", "CD4+ T cells")]
RNA.KEGG.sig.compare <- xEnrichCompare(RNA.sig.KEGG, FDR.cutoff = 0.01, 
                                       signature = FALSE, 
                                       bar.label = FALSE, 
                                       displayBy = "adjp")
RNA.KEGG.sig.compare <- 
    RNA.KEGG.sig.compare + theme(axis.text.y=element_text(size=12,color="black")) +
    scale_fill_manual(values = c(Okabe_Ito[3], Okabe_Ito[5], Okabe_Ito[6])) 
RNA.KEGG.sig.compare

dir.create("figures", showWarnings = FALSE)
ggsave("RNA.sig.KEGG.pdf", RNA.KEGG.sig.compare, width = 12, height = 14, useDingbats = FALSE, path = "./figures/")
ggsave("RNA.sig.KEGG.png", RNA.KEGG.sig.compare, width = 12, height = 14, device = "png", path = "./figures/")



###Repeat XGR analysis using cell-type specific backgrounds.

dir.create("GOBP.specific_background")
RNA.CD4.sig.GOBP.specific_background <- XGR_function(RNA.CD4.sig, RNA.CD4.background, 
                                                     ontology = "GOBP")
RNA.CD4.sig.GOBP.specific_background.table <- xEnrichViewer(RNA.CD4.sig.GOBP.specific_background, details = TRUE, 
                                        top_num = 50)
write.table(RNA.CD4.sig.GOBP.specific_background.table, "./GOBP.specific_background/RNA.CD4.sig.GOBP.specific_background.txt", sep = "\t", quote = FALSE)

RNA.CD8.sig.GOBP.specific_background <- XGR_function(RNA.CD8.sig, RNA.CD8.background, 
                                                     ontology = "GOBP")
RNA.CD8.sig.GOBP.specific_background.table <- xEnrichViewer(RNA.CD8.sig.GOBP.specific_background, details = TRUE, 
                                        top_num = 50)
write.table(RNA.CD8.sig.GOBP.specific_background.table, "./GOBP.specific_background/RNA.CD8.sig.GOBP.specific_background.txt", sep = "\t", quote = FALSE)

RNA.CD14.sig.GOBP.specific_background <- XGR_function(RNA.CD14.sig, RNA.background, 
                                                      ontology = "GOBP")
RNA.CD14.sig.GOBP.specific_background.table <- xEnrichViewer(RNA.CD14.sig.GOBP.specific_background, details = TRUE, 
                                         top_num = 50)
write.table(RNA.CD14.sig.GOBP.specific_background.table, "./GOBP.specific_background/RNA.CD14.sig.GOBP.specific_background.txt", sep = "\t", quote = FALSE)

#Compare 3 cell types
RNA.sig.GOBP.specific_background <- list(RNA.CD4.sig.GOBP.specific_background, RNA.CD8.sig.GOBP.specific_background, RNA.CD14.sig.GOBP.specific_background)
names(RNA.sig.GOBP.specific_background) <- c("CD4", "CD8", "CD14")
#Set the order
RNA.sig.GOBP.specific_background <- RNA.sig.GOBP.specific_background[c("CD14", "CD8", "CD4")]
RNA.GOBP.specific_background.sig.compare <- xEnrichCompare(RNA.sig.GOBP.specific_background, 
                                                           FDR.cutoff = 0.05, facet = TRUE)
RNA.GOBP.specific_background.sig.compare
Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000") 
RNA.GOBP.specific_background.sig.compare <- 
    RNA.GOBP.specific_background.sig.compare + 
    theme(axis.text.y=element_text(size=14,color="black")) +
    scale_fill_manual(values = c(Okabe_Ito[3], Okabe_Ito[5], Okabe_Ito[6])) 
RNA.GOBP.specific_background.sig.compare

dir.create("figures", showWarnings = FALSE)
ggsave("RNA.sig.GOBP.specific_background.pdf", RNA.GOBP.specific_background.sig.compare, width = 15, height = 15, useDingbats = FALSE, path = "./figures/")
ggsave("RNA.sig.GOBP.specific_background.png", RNA.GOBP.specific_background.sig.compare, width = 25, height = 15, device = "png", path = "./figures/")

#Make heatmap of results 
#Read in results tables so don't have to re-run all above code
RNA.CD4.sig.GOBP.specific_background.table <- read.table ("./GOBP.specific_background/RNA.CD4.sig.GOBP.specific_background.txt", 
                                      header = TRUE, sep = "\t")
RNA.CD8.sig.GOBP.specific_background.table <- read.table ("./GOBP.specific_background/RNA.CD8.sig.GOBP.specific_background.txt", 
                                      header = TRUE, sep = "\t")
RNA.CD14.GOBP.specific_background.sig.table <- read.table ("./GOBP.specific_background/RNA.CD14.sig.GOBP.specific_background.txt", 
                                       header = TRUE, sep = "\t")

#First collect all the bits of data into a df then convert to a matrix
RNA.CD4.sig.GOBP.specific_background.padj <- data.frame(RNA.CD4.sig.GOBP.specific_background.table$name, 
                                    RNA.CD4.sig.GOBP.specific_background.table$adjp)
colnames(RNA.CD4.sig.GOBP.specific_background.padj) <- c("name", "CD4.padj")
RNA.CD4.sig.GOBP.specific_background.padj <- RNA.CD4.sig.GOBP.specific_background.padj %>% filter(CD4.padj < 0.01)

RNA.CD8.sig.GOBP.specific_background.padj <- data.frame(RNA.CD8.sig.GOBP.specific_background.table$name, 
                                    RNA.CD8.sig.GOBP.specific_background.table$adjp)
colnames(RNA.CD8.sig.GOBP.specific_background.padj) <- c("name", "CD8.padj")
RNA.CD8.sig.GOBP.specific_background.padj <- RNA.CD8.sig.GOBP.specific_background.padj %>% filter(CD8.padj < 0.01)

RNA.CD14.sig.GOBP.specific_background.padj <- data.frame(RNA.CD14.sig.GOBP.specific_background.table$name, 
                                     RNA.CD14.sig.GOBP.specific_background.table$adjp)
colnames(RNA.CD14.sig.GOBP.specific_background.padj) <- c("name", "CD14.padj")
RNA.CD14.sig.GOBP.specific_background.padj <- RNA.CD14.sig.GOBP.specific_background.padj %>% filter(CD14.padj < 0.01)

RNA.sig.GOBP.specific_background.padj <- RNA.CD4.sig.GOBP.specific_background.padj %>%
    full_join (RNA.CD8.sig.GOBP.specific_background.padj) %>%
    full_join (RNA.CD14.sig.GOBP.specific_background.padj) %>%
    #Add columns for -log10 padj
    mutate(CD4 = -log10(CD4.padj)) %>%
    mutate(CD8 = -log10(CD8.padj)) %>%
    mutate(CD14 = -log10(CD14.padj)) %>%
    dplyr::select (c(1, 5, 6, 7))

head(RNA.sig.GOBP.specific_background.padj)
#Table of -log10(padj) --> make heatmap in GraphPad Prism!
write.table(RNA.sig.GOBP.specific_background.padj, 
            "./GOBP.specific_background/RNA.sig.GOBP.specific_background.padj.txt", 
            quote = FALSE, 
            sep = "\t", 
            row.names = FALSE, 
            na = "")

#Do XGR using Reactome with specific background
dir.create("Reactome.specific_background")
RNA.CD4.sig.Reactome.specific_background <- XGR_function(RNA.CD4.sig, RNA.CD4.background, 
                                                     ontology = "MsigdbC2REACTOME")
RNA.CD4.sig.Reactome.specific_background.table <- xEnrichViewer(RNA.CD4.sig.Reactome.specific_background, details = TRUE, 
                                                            top_num = 50)
write.table(RNA.CD4.sig.Reactome.specific_background.table, "./Reactome.specific_background/RNA.CD4.sig.Reactome.specific_background.txt", sep = "\t", quote = FALSE)

RNA.CD8.sig.Reactome.specific_background <- XGR_function(RNA.CD8.sig, RNA.CD8.background, 
                                                     ontology = "MsigdbC2REACTOME")
RNA.CD8.sig.Reactome.specific_background.table <- xEnrichViewer(RNA.CD8.sig.Reactome.specific_background, details = TRUE, 
                                                            top_num = 50)
write.table(RNA.CD8.sig.Reactome.specific_background.table, "./Reactome.specific_background/RNA.CD8.sig.Reactome.specific_background.txt", sep = "\t", quote = FALSE)

RNA.CD14.sig.Reactome.specific_background <- XGR_function(RNA.CD14.sig, RNA.background, 
                                                      ontology = "MsigdbC2REACTOME")
RNA.CD14.sig.Reactome.specific_background.table <- xEnrichViewer(RNA.CD14.sig.Reactome.specific_background, details = TRUE, 
                                                             top_num = 50)
write.table(RNA.CD14.sig.Reactome.specific_background.table, "./Reactome.specific_background/RNA.CD14.sig.Reactome.specific_background.txt", sep = "\t", quote = FALSE)

#Compare 3 cell types
RNA.sig.Reactome.specific_background <- list(RNA.CD4.sig.Reactome.specific_background, RNA.CD8.sig.Reactome.specific_background, RNA.CD14.sig.Reactome.specific_background)
names(RNA.sig.Reactome.specific_background) <- c("CD4", "CD8", "CD14")
#Set the order
RNA.sig.Reactome.specific_background <- RNA.sig.Reactome.specific_background[c("CD14", "CD8", "CD4")]
RNA.Reactome.specific_background.sig.compare <- xEnrichCompare(RNA.sig.Reactome.specific_background, FDR.cutoff = 0.05, facet = TRUE)
RNA.Reactome.specific_background.sig.compare
Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000") 
RNA.Reactome.specific_background.sig.compare <- 
    RNA.Reactome.specific_background.sig.compare + theme(axis.text.y=element_text(size=14,color="black")) +
    scale_fill_manual(values = c(Okabe_Ito[3], Okabe_Ito[5], Okabe_Ito[6])) 
RNA.Reactome.specific_background.sig.compare

dir.create("figures", showWarnings = FALSE)
ggsave("RNA.sig.Reactome.specific_background.pdf", RNA.Reactome.specific_background.sig.compare, width = 15, height = 15, useDingbats = FALSE, path = "./figures/")
ggsave("RNA.sig.Reactome.specific_background.png", RNA.Reactome.specific_background.sig.compare, width = 25, height = 15, device = "png", path = "./figures/")

#Make heatmap of results 
#Read in results tables so don't have to re-run all above code
RNA.CD4.sig.Reactome.specific_background.table <- read.table ("./Reactome.specific_background/RNA.CD4.sig.Reactome.specific_background.txt", 
                                                          header = TRUE, sep = "\t")
RNA.CD8.sig.Reactome.specific_background.table <- read.table ("./Reactome.specific_background/RNA.CD8.sig.Reactome.specific_background.txt", 
                                                          header = TRUE, sep = "\t")
RNA.CD14.Reactome.specific_background.sig.table <- read.table ("./Reactome.specific_background/RNA.CD14.sig.Reactome.specific_background.txt", 
                                                           header = TRUE, sep = "\t")

#First collect all the bits of data into a df then convert to a matrix
RNA.CD4.sig.Reactome.specific_background.padj <- data.frame(RNA.CD4.sig.Reactome.specific_background.table$name, 
                                                        RNA.CD4.sig.Reactome.specific_background.table$adjp)
colnames(RNA.CD4.sig.Reactome.specific_background.padj) <- c("name", "CD4.padj")
RNA.CD4.sig.Reactome.specific_background.padj <- RNA.CD4.sig.Reactome.specific_background.padj %>% filter(CD4.padj < 0.01)

RNA.CD8.sig.Reactome.specific_background.padj <- data.frame(RNA.CD8.sig.Reactome.specific_background.table$name, 
                                                        RNA.CD8.sig.Reactome.specific_background.table$adjp)
colnames(RNA.CD8.sig.Reactome.specific_background.padj) <- c("name", "CD8.padj")
RNA.CD8.sig.Reactome.specific_background.padj <- RNA.CD8.sig.Reactome.specific_background.padj %>% filter(CD8.padj < 0.01)

RNA.CD14.sig.Reactome.specific_background.padj <- data.frame(RNA.CD14.sig.Reactome.specific_background.table$name, 
                                                         RNA.CD14.sig.Reactome.specific_background.table$adjp)
colnames(RNA.CD14.sig.Reactome.specific_background.padj) <- c("name", "CD14.padj")
RNA.CD14.sig.Reactome.specific_background.padj <- RNA.CD14.sig.Reactome.specific_background.padj %>% filter(CD14.padj < 0.01)

RNA.sig.Reactome.specific_background.padj <- RNA.CD4.sig.Reactome.specific_background.padj %>%
    full_join (RNA.CD8.sig.Reactome.specific_background.padj) %>%
    full_join (RNA.CD14.sig.Reactome.specific_background.padj) %>%
    #Add columns for -log10 padj
    mutate(CD4 = -log10(CD4.padj)) %>%
    mutate(CD8 = -log10(CD8.padj)) %>%
    mutate(CD14 = -log10(CD14.padj)) %>%
    dplyr::select (c(1, 5, 6, 7))

head(RNA.sig.Reactome.specific_background.padj)
#Table of -log10(padj) --> make heatmap in GraphPad Prism!
write.table(RNA.sig.Reactome.specific_background.padj, 
            "./Reactome.specific_background/RNA.sig.Reactome.specific_background.padj.txt", 
            quote = FALSE, 
            sep = "\t", 
            row.names = FALSE, 
            na = "")






###Are any DE genes near GWAS loci?###

#Create a function called get_gene_locations
# get locations of genes from members_Overlap using BioMart
# create a df in bed format of the output
library(biomaRt)
mart <- useEnsembl("genes", "hsapiens_gene_ensembl", GRCh = "37")

get_gene_locations <- function(df){
    df.loc <- getBM(attributes=c('chromosome_name', 
                                  'start_position',
                                  'end_position', 
                                  'hgnc_symbol'), 
                     filters = 'hgnc_symbol', 
                     values = df$name, 
                     mart = mart)
    colnames(df.loc) <- c("chr", "start", "end", "gene")
    #remove columns where chr is not just a digit
    df.loc <- subset (df.loc, grepl("^[[:digit:]]", df.loc$chr))
    #Add chr to start of the chr column
    df.loc <- df.loc %>% mutate (chr = paste0("chr", chr))
    print(head(df.loc))
    return(df.loc)
}


#create a list of inputs then run the function over the list

RNA.sig.genes.bed <- lapply (RNA.sig, function (x) get_gene_locations(x))
names(RNA.sig.genes.bed) <- c("RNA.CD4.sig", "RNA.CD8.sig", "RNA.CD14.sig")

dir.create("gene_locations")
for (i in 1:length(RNA.sig.genes.bed)){
    write.table (RNA.sig.genes.bed[[i]], 
                 paste0("./gene_locations/",names(RNA.sig.genes.bed[i]),".bed"),
                 col.names = FALSE, 
                 row.names = FALSE,
                 quote = FALSE,
                 sep = "\t")
}

#Save all the eTerms for later use
save(RNA.CD4.sig.GOBP, RNA.CD8.sig.GOBP, RNA.CD14.sig.GOBP, 
     RNA.CD4.sig.Reactome, RNA.CD8.sig.Reactome, RNA.CD14.sig.Reactome, file = "./eTerms/eTerms.RData" )

save(RNA.CD4.sig.GOBP.specific_background, RNA.CD8.sig.GOBP.specific_background, RNA.CD14.sig.GOBP.specific_background, 
     RNA.CD4.sig.Reactome.specific_background, RNA.CD8.sig.Reactome.specific_background, RNA.CD14.sig.Reactome.specific_background, 
     file = "./eTerms/eTerms_specific_background.RData" )
