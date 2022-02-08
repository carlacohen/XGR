#Analysis on RNA-seq cohorts 2-5

setwd("~/R/XGR/RNAseq_cohort2-5")
#Aim 
# Run XGR on DE genes

library(dplyr)
library(stringr)
library(tidyr)
library(XGR)
library(ggplot2)

#Define colours
Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000") 


#Read in RNA data (output from DEseq2)
RNA.CD4 <- read.table("Inputs/RNA.CD4.txt", header = TRUE, sep = "\t") 
RNA.CD8 <- read.table("Inputs/RNA.CD8.txt", header = TRUE, sep = "\t") 
RNA.CD14 <- read.table("Inputs/RNA.CD14.txt", header = TRUE, sep = "\t") 

#Generate df of gene with padj score for each cell type

RNA.CD4.genes <- RNA.CD4 %>%
    select(name, padj)

RNA.CD8.genes <- RNA.CD8 %>%
    select(name, padj)

RNA.CD14.genes <- RNA.CD14 %>%
    select(name, padj)

#Performe XGR network analysis

dir.create("Networks")
RData.location <- "http://galahad.well.ox.ac.uk/bigdata"

XGR_network <- function(input){
    input_network <- xSubneterGenes(data = input,
                                    network = "STRING_high",
                                    subnet.significance=0.01,
                                    subnet.size=50,
                                    RData.location = RData.location)
    return(input_network)
}

RNA.CD4.network <- XGR_network(RNA.CD4.genes)
pattern.CD4 <- -log10(as.numeric(V(RNA.CD4.network)$significance))
RNA.CD4.network.plot <- xVisNet(g=RNA.CD4.network, pattern=pattern.CD4, vertex.shape="sphere", vertex.label.font=2, newpage=F)


RNA.CD8.network <- XGR_network(RNA.CD8.genes)
pattern.CD8 <- -log10(as.numeric(V(RNA.CD8.network)$significance))
RNA.CD8.network.plot <- xVisNet(g=RNA.CD8.network, pattern=pattern.CD8, vertex.shape="sphere", vertex.label.font=2, newpage=F)
RNA.CD8.network.plot

RNA.CD14.network <- XGR_network(RNA.CD14.genes)
pattern.CD14 <- -log10(as.numeric(V(RNA.CD14.network)$significance))
RNA.CD14.network.plot <- xVisNet(g=RNA.CD14.network, pattern=pattern.CD14, vertex.shape="sphere", vertex.label.font=2, newpage=F)
