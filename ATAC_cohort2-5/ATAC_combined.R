#Aim 
#generate lists of genes assigned to ATAC peaks via
#1 proximal gene
#2 PCHIC
#3 overlap with TSS
# Then run XGR on the combined set of genes

library(dplyr)
library(stringr)
library(tidyr)

Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000") 


#Read in ATAC data (output from DEseq2)
ATAC.CD4 <- read.table("Inputs/ATAC.CD4.txt", header = TRUE, row.names = NULL, sep = "\t") 
ATAC.CD8 <- read.table("Inputs/ATAC.CD8.txt", header = TRUE, row.names = NULL, sep = "\t") 
ATAC.CD14 <- read.table("Inputs/ATAC.CD14.txt", header = TRUE, row.names = NULL, sep = "\t") 

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


#How many peaks have FC > 1.5 and padj < 0.05?
sig_peaks <- function (df) {
    df2 <- df %>% 
        filter(log2FoldChange > 0.585 | log2FoldChange < -0.585) %>%
        filter(padj<0.05) 
    return (df2)
}

#Run function on all 3 cell type
ATAC.CD4.sig <- sig_peaks(ATAC.CD4)
ATAC.CD8.sig <- sig_peaks(ATAC.CD8)
ATAC.CD14.sig <- sig_peaks(ATAC.CD14)
dim(ATAC.CD4.sig) # 0 peaks
dim (ATAC.CD8.sig) #0 peaks
dim (ATAC.CD14.sig) #1 peaks

#Need to use top 200 peaks for pathway analysis.

#Create a function to select top 200 peaks with FC > 1.5 and sorted by p value
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
dim (ATAC.CD4.Top200) #97 peaks
dim (ATAC.CD8.Top200) #200 peaks
dim (ATAC.CD14.Top200) #84 peaks

#Create a bed file of Top 200 significant regions

#4th column is gene name
make_bed <- function (df) {
    df2 <- df %>% 
        dplyr::select (chr, start, end, ProxGene) %>%
        mutate_all(na_if,"") %>%
        #remove rows where col 3 == col 2
        filter (start != end)
    return (df2)
}

ATAC.CD4.Top200.bed <- make_bed(ATAC.CD4.Top200)
ATAC.CD8.Top200.bed <- make_bed(ATAC.CD8.Top200)
ATAC.CD14.Top200.bed <- make_bed(ATAC.CD14.Top200)

#4th column is ones
make_bed_2 <- function (df) {
    df2 <- df %>% 
        dplyr::select (chr, start, end) %>%
        mutate_all(na_if,"") %>%
        #remove rows where col 3 == col 2
        filter (start != end) %>%
        #add a column of ones
        mutate("1" =1)
    return (df2)
}

ATAC.CD4.Top200.bed <- make_bed_2(ATAC.CD4.Top200)
ATAC.CD8.Top200.bed <- make_bed_2(ATAC.CD8.Top200)
ATAC.CD14.Top200.bed <- make_bed_2(ATAC.CD14.Top200)

#save bed files

dir.create("bed")
write.table (ATAC.CD4.Top200.bed, "./bed/ATAC.CD4.Top200.bed", 
             col.names = FALSE,
             row.names = FALSE,
             quote = FALSE,
             na = "NA",
             sep = "\t")
write.table (ATAC.CD8.Top200.bed, "./bed/ATAC.CD8.Top200.bed", 
             col.names = FALSE,
             row.names = FALSE,
             quote = FALSE,
             na = "NA",
             sep = "\t")
write.table (ATAC.CD14.Top200.bed, "./bed/ATAC.CD14.Top200.bed", 
             col.names = FALSE,
             row.names = FALSE,
             quote = FALSE,
             na = "NA",
             sep = "\t")

#on galahad at http://galahad.well.ox.ac.uk/UCSC/Carla_WashU/Differential_peaks/



#Get list of genes assigned by proximity
ATAC.CD4.Top200.proximity <- ATAC.CD4.Top200 %>% 
    select (ProxGene) %>% 
    rename("gene" = "ProxGene")
ATAC.CD8.Top200.proximity <- ATAC.CD8.Top200 %>% 
    select (ProxGene) %>% 
    rename("gene" = "ProxGene")
ATAC.CD14.Top200.proximity <- ATAC.CD14.Top200 %>% 
    select (ProxGene) %>% 
    rename("gene" = "ProxGene")



#Background lists of genes assigned by proximity
ATAC.CD4.background.proximity <- ATAC.CD4 %>% select (ProxGene) %>% rename("gene" = "ProxGene")
ATAC.CD8.background.proximity <- ATAC.CD8 %>% select (ProxGene) %>% rename("gene" = "ProxGene")
ATAC.CD14.background.proximity <- ATAC.CD14 %>% select (ProxGene) %>% rename("gene" = "ProxGene")

#Make a concatenated background set using rbind
ATAC.background.proximity <- rbind (ATAC.CD4.background.proximity, 
                                ATAC.CD8.background.proximity, 
                                ATAC.CD14.background.proximity)


#Read in data of peaks linked to genes via PCHIC
ATAC.CD4.PCHIC <- read.table("Inputs/PCHIC/ATAC_ML_CD4_PCHIC.txt", header = FALSE, sep = "\t")
colnames (ATAC.CD4.PCHIC) <- c("chr", "start", "end", "1", "chr_loop", "start_loop", "end_loop", "interaction", "gene", "distance")
ATAC.CD8.PCHIC <- read.table("Inputs/PCHIC/ATAC_ML_CD8_PCHIC.txt", header = FALSE, sep = "\t")
colnames (ATAC.CD8.PCHIC) <- c("chr", "start", "end", "1", "chr_loop", "start_loop", "end_loop", "interaction", "gene", "distance")
ATAC.CD14.PCHIC <- read.table("Inputs/PCHIC/ATAC_ML_CD14_PCHIC.txt", header = FALSE, sep = "\t")
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
dim(ATAC.CD4.Top200.PCHIC) #37 genes
dim(ATAC.CD8.Top200.PCHIC) #72 genes
dim(ATAC.CD14.Top200.PCHIC) #23 genes

#Generate background set for PCHIC linked genes
extract.ATAC.background.PCHIC <- function (ATAC.ML.PCHIC){
    #Filter out rows where gene is "."
    ATAC.ML.PCHIC.genes <- ATAC.ML.PCHIC %>% filter (gene !=".")
    #split genes separated by ";" into a single list
    list <- strsplit(ATAC.ML.PCHIC.genes$gene, ";")
    # use unlist to put each gene on a new row
    final_list <- as.data.frame(unlist(list))
    #show only unique genes
    unique_list <- unique(final_list)
    #rename columns to "gene"
    colnames(unique_list) <- c("gene")
    return(unique_list)
}

ATAC.CD4.background.PCHIC <- extract.ATAC.background.PCHIC(ATAC.CD4.PCHIC)
ATAC.CD8.background.PCHIC <- extract.ATAC.background.PCHIC(ATAC.CD8.PCHIC)
ATAC.CD14.background.PCHIC <- extract.ATAC.background.PCHIC(ATAC.CD14.PCHIC)
#Make a concatenated background set using rbind
ATAC.background.PCHIC <- rbind (ATAC.CD4.background.PCHIC, 
                                ATAC.CD8.background.PCHIC, 
                                ATAC.CD14.background.PCHIC)

#Generate concatenated lists of genes from both methods
ATAC.CD4.genes <- rbind (ATAC.CD4.Top200.proximity, ATAC.CD4.Top200.PCHIC)
ATAC.CD8.genes <- rbind (ATAC.CD8.Top200.proximity, ATAC.CD8.Top200.PCHIC)
ATAC.CD14.genes <- rbind (ATAC.CD14.Top200.proximity, ATAC.CD14.Top200.PCHIC)
#Convert these lists to vectors
ATAC.CD4.genes <- as.vector(ATAC.CD4.genes$gene)
ATAC.CD8.genes <- as.vector(ATAC.CD8.genes$gene)
ATAC.CD14.genes <- as.vector(ATAC.CD14.genes$gene)
length(ATAC.CD4.genes) #134 genes
length(ATAC.CD8.genes) #272 genes
length(ATAC.CD14.genes) #107
#So I can use the default parameters for XGR for all 3 cell types

#Generate concatenated background for all 3 cell types
ATAC.background <- rbind(ATAC.background.proximity, ATAC.background.PCHIC)
ATAC.background <- as.vector(ATAC.background$gene)
class(ATAC.background)

library(XGR)
dir.create("GOBP")
dir.create("Reactome")
RData.location <- "http://galahad.well.ox.ac.uk/bigdata"

XGR_function <- function(genelist, background, ontology) {
    #Run xEnricher function 
    #NB settings: for <50 genes use size.range = c(5,2000) min.overlap = 3
    #For larger gene sets than 50 use size.range = c(10,2000) min.overlap = 5 (which is the default setting)
    enriched_genes <- xEnricherGenes(data = genelist, 
                                     background = background, 
                                     test = "hypergeo", 
                                     ontology = ontology, 
                                     RData.location = RData.location)
    enriched_genes_concise <- xEnrichConciser(enriched_genes)
    return (enriched_genes_concise)
}

ATAC.CD4.GOBP <- XGR_function(ATAC.CD4.genes, ATAC.background, ontology = "GOBP")
ATAC.CD4.GOBP.table <- xEnrichViewer(ATAC.CD4.GOBP, details = TRUE, 
                                    top_num = 50)
write.table(ATAC.CD4.GOBP.table, "./GOBP/ATAC.CD4.GOBP.txt", sep = "\t", quote = FALSE)

ATAC.CD8.GOBP <- XGR_function(ATAC.CD8.genes, ATAC.background, ontology = "GOBP")
ATAC.CD8.GOBP.table <- xEnrichViewer(ATAC.CD8.GOBP, details = TRUE, 
                                     top_num = 50)
write.table(ATAC.CD8.GOBP.table, "./GOBP/ATAC.CD8.GOBP.txt", sep = "\t", quote = FALSE)

ATAC.CD14.GOBP <- XGR_function(ATAC.CD14.genes, ATAC.background, ontology = "GOBP")
ATAC.CD14.GOBP.table <- xEnrichViewer(ATAC.CD14.GOBP, details = TRUE, 
                                     top_num = 50)
write.table(ATAC.CD14.GOBP.table, "./GOBP/ATAC.CD14.GOBP.txt", sep = "\t", quote = FALSE)

#make a list of things to compare
ATAC.GOBP <- list(ATAC.CD4.GOBP, ATAC.CD8.GOBP, ATAC.CD14.GOBP)
names(ATAC.GOBP) <- c("CD4 T cells", "CD8 T cells", "CD14 monocytes")

#Set the order
ATAC.GOBP <- ATAC.GOBP[c("CD14 monocytes", "CD8 T cells", "CD4 T cells")]

ATAC.GOBP.compare <- xEnrichCompare(ATAC.GOBP, FDR.cutoff = 0.05, facet = TRUE, signature = FALSE,
                                    bar.label.size = 4)
ATAC.GOBP.compare <- 
    ATAC.GOBP.compare + theme(axis.text.y=element_text(size=20,color="black")) +
    scale_fill_manual(values = c(Okabe_Ito[3], Okabe_Ito[5], Okabe_Ito[6])) 
ATAC.GOBP.compare

dir.create("figures", showWarnings = FALSE)
ggsave("ATAC.GOBP.pdf", ATAC.GOBP.compare, width = 15, height = 5, useDingbats = FALSE, path = "./figures/")
ggsave("ATAC.GOBP.png", ATAC.GOBP.compare, width = 15, height = 5, device = "png", path = "./figures/")




#Now do for Reactome

ATAC.CD4.Reactome <- XGR_function(ATAC.CD4.genes, ATAC.background, ontology = "MsigdbC2REACTOME")

ATAC.CD4.Reactome.table <- xEnrichViewer(ATAC.CD4.Reactome, details = TRUE, 
                                     top_num = 50)

write.table(ATAC.CD4.Reactome.table, "./Reactome/ATAC.CD4.Reactome.txt", sep = "\t", quote = FALSE)

ATAC.CD8.Reactome <- XGR_function(ATAC.CD8.genes, ATAC.background, ontology = "MsigdbC2REACTOME")
ATAC.CD8.Reactome.table <- xEnrichViewer(ATAC.CD8.Reactome, details = TRUE, 
                                     top_num = 50)
write.table(ATAC.CD8.Reactome.table, "./Reactome/ATAC.CD8.Reactome.txt", sep = "\t", quote = FALSE)

ATAC.CD14.Reactome <- XGR_function(ATAC.CD14.genes, ATAC.background, ontology = "MsigdbC2REACTOME")
ATAC.CD14.Reactome.table <- xEnrichViewer(ATAC.CD14.Reactome, details = TRUE, 
                                      top_num = 50)
write.table(ATAC.CD14.Reactome.table, "./Reactome/ATAC.CD14.Reactome.txt", sep = "\t", quote = FALSE)

#Perform comparison
ATAC.Reactome <- list(ATAC.CD4.Reactome, ATAC.CD8.Reactome, ATAC.CD14.Reactome)
names(ATAC.Reactome) <- c("CD4+ T cells", "CD8+ T cells", "CD14+ monocytes")
#Set the order
ATAC.Reactome <- ATAC.Reactome[c("CD14+ monocytes", "CD8+ T cells", "CD4+ T cells")]

#Perform the comparison
ATAC.Reactome.compare <- xEnrichCompare(ATAC.Reactome, FDR.cutoff = 0.05, 
                                           signature = FALSE,
                                           bar.label.size = 4)

ATAC.Reactome.compare <- 
    ATAC.Reactome.compare + theme(axis.text.y=element_text(size=24,color="black")) +
    scale_fill_manual(values = c(Okabe_Ito[3], Okabe_Ito[5], Okabe_Ito[6])) 
ATAC.Reactome.compare

dir.create("figures", showWarnings = FALSE)
ggsave("ATAC.Reactome.pdf", ATAC.Reactome.compare, width = 15, height = 10, useDingbats = FALSE, path = "./figures/")
ggsave("ATAC.Reactome.png", ATAC.Reactome.compare, width = 15, height = 10, device = "png", path = "./figures/")


#Repeat XGR analysis with less stringent parameters for Reactome
XGR_function_2 <- function(genelist, background, ontology) {
    #Run xEnricher function 
    #NB settings: for <50 genes use size.range = c(5,2000) min.overlap = 3
    #For larger gene sets than 50 use size.range = c(10,2000) min.overlap = 5 (which is the default setting)
    enriched_genes <- xEnricherGenes(data = genelist, 
                                     background = background, 
                                     test = "hypergeo", 
                                     ontology = ontology, 
                                     RData.location = RData.location, 
                                     size.range = c(5,2000), 
                                     min.overlap = 3)
    enriched_genes_concise <- xEnrichConciser(enriched_genes)
    return (enriched_genes_concise)
}

ATAC.CD4.Reactome_2 <- XGR_function_2(ATAC.CD4.genes, ATAC.background, ontology = "MsigdbC2REACTOME")

ATAC.CD4.Reactome.table_2 <- xEnrichViewer(ATAC.CD4.Reactome_2, details = TRUE, 
                                         top_num = 50)

write.table(ATAC.CD4.Reactome.table_2, "./Reactome/ATAC.CD4.Reactome_2.txt", sep = "\t", quote = FALSE)

ATAC.CD8.Reactome_2 <- XGR_function_2(ATAC.CD8.genes, ATAC.background, ontology = "MsigdbC2REACTOME")

ATAC.CD8.Reactome.table_2 <- xEnrichViewer(ATAC.CD8.Reactome_2, details = TRUE, 
                                           top_num = 50)

write.table(ATAC.CD8.Reactome.table_2, "./Reactome/ATAC.CD4.Reactome_2.txt", sep = "\t", quote = FALSE)

ATAC.CD14.Reactome_2 <- XGR_function_2(ATAC.CD14.genes, ATAC.background, ontology = "MsigdbC2REACTOME")

ATAC.CD14.Reactome.table_2 <- xEnrichViewer(ATAC.CD14.Reactome_2, details = TRUE, 
                                           top_num = 50)

write.table(ATAC.CD4.Reactome.table_2, "./Reactome/ATAC.CD4.Reactome_2.txt", sep = "\t", quote = FALSE)
write.table(ATAC.CD8.Reactome.table_2, "./Reactome/ATAC.CD8.Reactome_2.txt", sep = "\t", quote = FALSE)
write.table(ATAC.CD14.Reactome.table_2, "./Reactome/ATAC.CD14.Reactome_2.txt", sep = "\t", quote = FALSE)


#Perform comparison
ATAC.Reactome_2 <- list(ATAC.CD4.Reactome_2, ATAC.CD8.Reactome_2, ATAC.CD14.Reactome_2)
names(ATAC.Reactome_2) <- c("CD4+ T cells", "CD8+ T cells", "CD14+ monocytes")
#Set the order
ATAC.Reactome_2 <- ATAC.Reactome_2[c("CD14+ monocytes", "CD8+ T cells", "CD4+ T cells")]

#Perform the comparison
ATAC.Reactome.compare_2 <- xEnrichCompare(ATAC.Reactome_2, FDR.cutoff = 0.05, 
                                        signature = FALSE,
                                        bar.label.size = 4)

ATAC.Reactome.compare_2 <- 
    ATAC.Reactome.compare_2 + theme(axis.text.y=element_text(size = 16, color="black")) +
    scale_fill_manual(values = c(Okabe_Ito[3], Okabe_Ito[5], Okabe_Ito[6])) 
ATAC.Reactome.compare_2

dir.create("figures", showWarnings = FALSE)
ggsave("ATAC.Reactome_2.pdf", ATAC.Reactome.compare_2, width = 15, height = 10, useDingbats = FALSE, path = "./figures/")
ggsave("ATAC.Reactome_2.png", ATAC.Reactome.compare_2, width = 15, height = 10, device = "png", path = "./figures/")


#Find examples from output 
# e.g. CCR4 is in top Reactome category for CD4
#Is CCR4 in the list of PCHIC genes or proximity genes?
ATAC.CD4.Top200.PCHIC %>% filter (gene == "CCR4")
#No
ATAC.CD4.Top200.proximity %>% filter (gene == "CCR4")
#Yes

#CCR7
#Is CCR7 in the list of PCHIC genes or proximity genes?
ATAC.CD4.Top200.PCHIC %>% filter (gene == "CCR7")
#No
ATAC.CD4.Top200.proximity %>% filter (gene == "CCR7")
#Yes


#GRK5
#Is GRK5 in the list of PCHIC genes or proximity genes?
ATAC.CD4.Top200.PCHIC %>% filter (gene == "GRK5")
#Yes
ATAC.CD4.Top200.proximity %>% filter (gene == "GRK5")
#Yes

#Find location of ATAC peak in PCHIC set & generate interaction track file of ATAC peaks
GRK5_loc <- ATAC.CD14.PCHIC %>% filter (grepl("GRK5", gene)) %>%
    select(chr, start, end, interaction)
GRK5_loc
dir.create("Examples")
write.table(GRK5_loc, "./Examples/ATAC_GRK5_interaction.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

#Find the p values of these regions
GRK5_info <- GRK5_loc %>% inner_join(ATAC.CD14) %>%
    arrange (padj)
GRK5_info
#One ATAC peak at chr10 121166080 121167442 is significant p<0.05 log2FC -0.17

#generate interaction file for text upload
#format chr1,713605,715737     chr1,720589,722848      2
library(stringr)
GRK5_loc_text <- ATAC.CD14.PCHIC %>% filter (grepl("GRK5", gene)) %>%
    select(chr, start, end, interaction)
GRK5_loc_text$ATAC_loc <- paste(GRK5_loc_text$chr, GRK5_loc_text$start, GRK5_loc_text$end, sep = ",")
GRK5_loc_text
gene_loc <- str_split(GRK5_loc_text$interaction, ",")
gene_loc
GRK5_loc_text$gene_loc <- lapply(gene_loc, function(x) x[1])
GRK5_loc_text$score <- lapply(gene_loc, function (x) x[2])
GRK5_loc_text
GRK5_loc_text$gene_loc <- unlist(GRK5_loc_text$gene_loc)
GRK5_loc_text$score <- unlist(GRK5_loc_text$score)
GRK5_loc_text <- GRK5_loc_text %>% select (ATAC_loc, gene_loc, score)
class(GRK5_loc_text)
write.table (GRK5_loc_text, "ATAC_GRK5_interaction_local.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

#Look in glutathione metabolism in CD14 as could be interesting pathway
#GSTA1, GSTM1, GSTP1

#GSTA1
#Is GSTA1 in the list of PCHIC genes or proximity genes?
ATAC.CD14.Top200.PCHIC %>% filter (gene == "GSTA1")
#No
ATAC.CD14.Top200.proximity %>% filter (gene == "GSTA1")
#Yes

#GSTM1
#Is GSTM1 in the list of PCHIC genes or proximity genes?
ATAC.CD14.Top200.PCHIC %>% filter (gene == "GSTM1")
#No
ATAC.CD14.Top200.proximity %>% filter (gene == "GSTM1")
#Yes

#GSTP1
#Is GSTP1 in the list of PCHIC genes or proximity genes?
ATAC.CD14.Top200.PCHIC %>% filter (gene == "GSTP1")
#No
ATAC.CD14.Top200.proximity %>% filter (gene == "GSTP1")
#Yes



#Do any top 200 genes also have significantly different gene expression?

#Import DEseq2 data from RNA-seq

RNA.CD4 <- read.table("C:/Users/ccohen.WHGC/Documents/R/XGR/RNAseq_cohort2-5/Inputs/RNA.CD4.txt",
                      header = TRUE, sep = "\t") 
RNA.CD8 <- read.table("C:/Users/ccohen.WHGC/Documents/R/XGR/RNAseq_cohort2-5/Inputs/RNA.CD8.txt",
                      header = TRUE, sep = "\t") 
RNA.CD14 <- read.table("C:/Users/ccohen.WHGC/Documents/R/XGR/RNAseq_cohort2-5/Inputs/RNA.CD14.txt",
                      header = TRUE, sep = "\t") 
#Intersect Top 200 ATAC peaks with RNAseq data and select RNA >1.5 fold change and padj<0.05
CD4.ATAC.RNA.sig <- ATAC.CD4.Top200 %>% inner_join(RNA.CD4, by = c("ProxGene" = "name")) %>%
    filter(log2FoldChange.y > 0.585 | log2FoldChange.y < -0.585) %>%
    filter(padj.y<0.05) 
CD8.ATAC.RNA.sig <- ATAC.CD8.Top200 %>% inner_join(RNA.CD8, by = c("ProxGene" = "name")) %>%
    filter(log2FoldChange.y > 0.585 | log2FoldChange.y < -0.585) %>%
    filter(padj.y<0.05) 
CD14.ATAC.RNA.sig <- ATAC.CD14.Top200 %>% inner_join(RNA.CD14, by = c("ProxGene" = "name")) %>%
    filter(log2FoldChange.y > 0.585 | log2FoldChange.y < -0.585) %>%
    filter(padj.y<0.05) 
#Yes there are a limited number 


#CD4
#HDGFL3 has biggest fold change 0.7350390
#CD8
#TNFSF14 interesting gene, FC 0.6835527
#CD14
#AHRR FC 0.61 (2 peaks)
#OLFM2 0.61



###Are any genes from enriched pathways near GWAS loci?### Not using
#Perhaps better at this stage to see if there are any Top200 DE peaks near SNPs?

#Create a function called get_gene_locations
# get locations of genes from members_Overlap using BioMart
# create a df in bed format of the output
library(biomaRt)
mart <- useEnsembl("genes", "hsapiens_gene_ensembl", GRCh = "37")

head(ATAC.CD8.GOBP.table[,13])
get_gene_locations <- function(df){
    df2 <- str_split(df[,13], ",")
    df2 <- unlist(df2)
    df2 <- unique(df2)
    df2 <- trimws(df2)
    print(head(df2))
    df2.loc <- getBM(attributes=c('chromosome_name', 
                                  'start_position',
                                  'end_position', 
                                  'hgnc_symbol'), 
                                  filters = 'hgnc_symbol', 
                                  values = df2, 
                                  mart = mart)
    colnames(df2.loc) <- c("chr", "start", "end", "gene")
    #remove columns where chr is not just a digit
    df2.loc <- subset (df2.loc, grepl("^[[:digit:]]", df2.loc$chr))
    #Add chr to start of the chr column
    df2.loc <- df2.loc %>% mutate (chr = paste0("chr", chr))
    print(head(df2.loc))
    return(df2.loc)
}

#create a list of inputs then run the function over the list
gene_members_list <- list(ATAC.CD14.GOBP.table, ATAC.CD4.GOBP.table, ATAC.CD8.GOBP.table,
                          ATAC.CD14.Reactome.table_2, ATAC.CD4.Reactome.table_2, ATAC.CD8.Reactome.table_2)
gene_members_list_bed <- lapply(gene_members_list, function (x) get_gene_locations(x))
names(gene_members_list_bed) <- c("ATAC.CD14.GOBP", "ATAC.CD4.GOBP", "ATAC.CD8.GOBP",
                          "ATAC.CD14.Reactome", "ATAC.CD4.Reactome", "ATAC.CD8.Reactome")
head(gene_members_list_bed[[1]])
head(gene_members_list_bed[[2]])
head(gene_members_list_bed[[3]])

dir.create("gene_locations")
for (i in 1:length(gene_members_list_bed)){
    write.table (gene_members_list_bed[[i]], 
                 paste0("./gene_locations/",names(gene_members_list_bed[i]),"_gene_members.bed"),
                 col.names = FALSE, 
                 row.names = FALSE,
                 quote = FALSE,
                 sep = "\t")
    }

#Now use bedtools to intersect with SNPs
#See readme file for details

#Read in the interaction results

ATAC.CD4.SNPs <- read.table("Top200.SNPs/ATAC.CD4.Top200.SNPs.txt", header = FALSE, row.names = NULL, sep = "\t") 
colnames (ATAC.CD4.SNPs) <- c("snp.chr", "snp.start", "snp.end", "snp", "peak.chr", "peak.start", "peak.end", "1")
ATAC.CD8.SNPs <- read.table("Top200.SNPs/ATAC.CD8.Top200.SNPs.txt", header = FALSE, row.names = NULL, sep = "\t") 
colnames (ATAC.CD8.SNPs) <- c("snp.chr", "snp.start", "snp.end", "snp", "peak.chr", "peak.start", "peak.end", "1")
ATAC.CD14.SNPs <- read.table("Top200.SNPs/ATAC.CD14.Top200.SNPs.txt", header = FALSE, row.names = NULL, sep = "\t") 
colnames (ATAC.CD14.SNPs) <- c("snp.chr", "snp.start", "snp.end", "snp", "peak.chr", "peak.start", "peak.end", "1")

#CD14 hits
#chr15 67294268 67295367 no proximal gene or PCHIC gene
#chr19 10040993 10042234 proxixal to OLFM2



###Are any regions near GWAS loci?### Another method
#Used bedtools to look for any ATAC peaks within 500kb of a lead SNP
#Results are in ./peak_vs_SNP

#Read in the data

ATAC_ML_CD4.SNP <- read.table("./peaks_vs_SNPs/ATAC_ML_CD4.SNP.txt", header = FALSE, sep = "\t")
colnames (ATAC_ML_CD4.SNP) <- c("chr.SNP", "start.SNP", "end.SNP", "SNP", 
                                   "chr", "start", "end", "1")
ATAC_ML_CD8.SNP <- read.table("./peaks_vs_SNPs/ATAC_ML_CD8.SNP.txt", header = FALSE, sep = "\t")
colnames (ATAC_ML_CD8.SNP) <- c("chr.SNP", "start.SNP", "end.SNP", "SNP", 
                                   "chr", "start", "end", "1")
ATAC_ML_CD14.SNP <- read.table("./peaks_vs_SNPs/ATAC_ML_CD14.SNP.txt", header = FALSE, sep = "\t")
colnames (ATAC_ML_CD14.SNP) <- c("chr.SNP", "start.SNP", "end.SNP", "SNP", 
                                    "chr", "start", "end", "1")

#Join with DEseq2 output, select peaks with FC > 1.5, p<0.05 and sort by distance to SNP

SNPs_peaks <- function (df, DEseq2) {
    df2 <- df %>% left_join(DEseq2) %>%
        filter (pvalue < 0.05) %>%
        filter(log2FoldChange > 0.585 | log2FoldChange < -0.585) %>%
        mutate (snp_dist = end.SNP - end) 
    #change snp dist to always be positive
    df2$snp.dist <- abs(df2$snp_dist)    
    df2 <-  df2 %>% arrange (snp.dist) %>%
        dplyr::select (SNP, chr, start, end, baseMean, log2FoldChange, pvalue, padj, ProxGene, snp.dist)
    return (df2)
}

ATAC.CD4.SNP <- SNPs_peaks(ATAC_ML_CD4.SNP, ATAC.CD4)
ATAC.CD8.SNP <- SNPs_peaks(ATAC_ML_CD8.SNP, ATAC.CD8)
ATAC.CD14.SNP <- SNPs_peaks(ATAC_ML_CD14.SNP, ATAC.CD14)
dim(ATAC.CD4.SNP)
dim(ATAC.CD8.SNP)
dim(ATAC.CD14.SNP)

#Very few examples...2 in CD14
#One is the OLFM2 peak found above
#chr15 67294268 67295367 at rs35874463 near SMAD3

#CD4 chr5 158788839 158789901 at IL12B



write.table (ATAC.CD4.SNP, "./peaks_vs_SNPs/ATAC.CD4.SNP.txt", 
             col.names = TRUE,
             row.names = FALSE,
             quote = FALSE,
             na = "NA",
             sep = "\t")

write.table (ATAC.CD8.SNP, "./peaks_vs_SNPs/ATAC.CD8.SNP.txt", 
             col.names = TRUE,
             row.names = FALSE,
             quote = FALSE,
             na = "NA",
             sep = "\t")
write.table (ATAC.CD14.SNP, "./peaks_vs_SNPs/ATAC.CD14.SNP.txt", 
             col.names = TRUE,
             row.names = FALSE,
             quote = FALSE,
             na = "NA",
             sep = "\t")

#What are the top peaks and do they have sig RNA expression also?
head(ATAC.CD14.sig)

#Top peak in CD14 chr2 225183411 225184916 is near WDFY1, MRPL4, SERPINE2, FAM124B, CUL3
head(RNA.CD14)
RNA.CD14 %>% filter (name == "WDFY1" | name == "MRPL4" | name == "SERPINE2"|
                     name == "FAM124B" | name == "CUL3" ) %>%
    filter (pvalue < 0.05) %>%
    filter (log2FoldChange < -0.585 | log2FoldChange > 0.585)
#No genes have sig DE expression

#Top 200 CD14
head (ATAC.CD14.Top200)
# chr1 111329467 111331046
#Near CD53 (and eQTL), KCNA2, KCNA3, LRIF1
RNA.CD14 %>% filter (name == "CD53" | name == "KCNA2" | name == "KCNA3"|
                         name == "LRIF1") #%>%
    #filter (pvalue < 0.05) %>%
    #filter (log2FoldChange < -0.585 | log2FoldChange > 0.585)
#No genes have sig DE expression

ATAC.CD14.Top200 %>% 
    filter (pvalue < 0.05) %>% 
    filter (log2FoldChange < -0.8 | log2FoldChange > 0.8) %>%
    arrange(log2FoldChange)

#chr20   1575674   1576905 near SIRPB1
#Mapping is odd, actually looks like a number of smaller peaks. Also intronic.

#chr10  55824074  55825074 near in intron of PCDH15
#chr3 154842218 154843735 intron of MME
#Join Top200 ATAC peaks with DE genes...did that already.
#Continue when genome browser is working!


