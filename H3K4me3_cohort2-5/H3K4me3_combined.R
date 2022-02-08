#Aim 
#generate lists of genes assigned to H3K4me3 peaks via
#1 proximal gene
#2 PCHIC
#3 overlap with TSS
# Then run XGR on the combined set of genes
setwd("~/R/XGR/H3K4me3_cohort2-5")

library(dplyr)
library(stringr)
library(tidyr)

Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000") 


#Read in H3K4me3 data (output from DEseq2)
H3K4me3.CD4 <- read.table("Inputs/H3K4me3.CD4.txt", header = TRUE, row.names = NULL, sep = "\t") 
H3K4me3.CD8 <- read.table("Inputs/H3K4me3.CD8.txt", header = TRUE, row.names = NULL, sep = "\t") 
H3K4me3.CD14 <- read.table("Inputs/H3K4me3.CD14.txt", header = TRUE, row.names = NULL, sep = "\t") 

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


#Create a function to select top 200 peaks with FC > 1.5 and sorted by p value
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

#Create a bed file of Top 200 significant regions

make_bed <- function (df) {
    df2 <- df %>% 
        dplyr::select (chr, start, end, ProxGene) %>%
        mutate_all(na_if,"") %>%
        #remove rows where col 3 == col 2
        filter (start != end)
    return (df2)
}

H3K4me3.CD4.Top200.bed <- make_bed(H3K4me3.CD4.Top200)
H3K4me3.CD8.Top200.bed <- make_bed(H3K4me3.CD8.Top200)
H3K4me3.CD14.Top200.bed <- make_bed(H3K4me3.CD14.Top200)

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

H3K4me3.CD4.Top200.bed <- make_bed_2(H3K4me3.CD4.Top200)
H3K4me3.CD8.Top200.bed <- make_bed_2(H3K4me3.CD8.Top200)
H3K4me3.CD14.Top200.bed <- make_bed_2(H3K4me3.CD14.Top200)

#save bed files

dir.create("bed")
write.table (H3K4me3.CD4.Top200.bed, "./bed/H3K4me3.CD4.Top200.bed", 
                 col.names = FALSE,
                 row.names = FALSE,
                 quote = FALSE,
                 na = "NA",
                 sep = "\t")
write.table (H3K4me3.CD8.Top200.bed, "./bed/H3K4me3.CD8.Top200.bed", 
             col.names = FALSE,
             row.names = FALSE,
             quote = FALSE,
             na = "NA",
             sep = "\t")
write.table (H3K4me3.CD14.Top200.bed, "./bed/H3K4me3.CD14.Top200.bed", 
             col.names = FALSE,
             row.names = FALSE,
             quote = FALSE,
             na = "NA",
             sep = "\t")




#Get list of genes assigned by proximity
H3K4me3.CD4.Top200.proximity <- H3K4me3.CD4.Top200 %>% 
    select (ProxGene) %>% 
    rename("gene" = "ProxGene")
H3K4me3.CD8.Top200.proximity <- H3K4me3.CD8.Top200 %>% 
    select (ProxGene) %>% 
    rename("gene" = "ProxGene")
H3K4me3.CD14.Top200.proximity <- H3K4me3.CD14.Top200 %>% 
    select (ProxGene) %>% 
    rename("gene" = "ProxGene")

#Create a function to select top significant peaks with padj < 0.05,  FC > 1.5 
sig.peaks <- function (df) {
    df2 <- df %>% 
        filter(log2FoldChange > 0.585 | log2FoldChange < -0.585) %>%
        filter(padj < 0.05) 
    return (df2)
}

#Run function on all 3 cell type
H3K4me3.CD4.sig <- sig.peaks(H3K4me3.CD4)
H3K4me3.CD8.sig <- sig.peaks(H3K4me3.CD8)
H3K4me3.CD14.sig <- sig.peaks(H3K4me3.CD14)

#How many significant peaks?
dim(H3K4me3.CD4.sig) # CD4 278
dim (H3K4me3.CD8.sig) # CD8 1 peak
dim (H3K4me3.CD14.sig) # CD14 2 peaks
#so continue using top 200 peaks


#Get list of genes assigned by proximity
H3K4me3.CD4.Top200.proximity <- H3K4me3.CD4.Top200 %>% 
    select (ProxGene) %>% 
    rename("gene" = "ProxGene")
H3K4me3.CD8.Top200.proximity <- H3K4me3.CD8.Top200 %>% 
    select (ProxGene) %>% 
    rename("gene" = "ProxGene")
H3K4me3.CD14.Top200.proximity <- H3K4me3.CD14.Top200 %>% 
    select (ProxGene) %>% 
    rename("gene" = "ProxGene")





#Background lists of genes assigned by proximity
H3K4me3.CD4.background.proximity <- H3K4me3.CD4 %>% select (ProxGene) %>% rename("gene" = "ProxGene")
H3K4me3.CD8.background.proximity <- H3K4me3.CD8 %>% select (ProxGene) %>% rename("gene" = "ProxGene")
H3K4me3.CD14.background.proximity <- H3K4me3.CD14 %>% select (ProxGene) %>% rename("gene" = "ProxGene")

#Make a concatenated background set using rbind
H3K4me3.background.proximity <- rbind (H3K4me3.CD4.background.proximity, 
                                    H3K4me3.CD8.background.proximity, 
                                    H3K4me3.CD14.background.proximity)


#Read in data of peaks linked to genes via PCHIC
H3K4me3.CD4.PCHIC <- read.table("Inputs/PCHIC/H3K4me3_ML_CD4_PCHIC.txt", header = FALSE, sep = "\t")
colnames (H3K4me3.CD4.PCHIC) <- c("chr", "start", "end", "1", "chr_loop", "start_loop", "end_loop", "interaction", "gene", "distance")
H3K4me3.CD8.PCHIC <- read.table("Inputs/PCHIC/H3K4me3_ML_CD8_PCHIC.txt", header = FALSE, sep = "\t")
colnames (H3K4me3.CD8.PCHIC) <- c("chr", "start", "end", "1", "chr_loop", "start_loop", "end_loop", "interaction", "gene", "distance")
H3K4me3.CD14.PCHIC <- read.table("Inputs/PCHIC/H3K4me3_ML_CD14_PCHIC.txt", header = FALSE, sep = "\t")
colnames (H3K4me3.CD14.PCHIC) <- c("chr", "start", "end","1", "chr_loop", "start_loop", "end_loop", "interaction", "gene", "distance")

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

#Generate background set for PCHIC linked genes
extract.H3K4me3.background.PCHIC <- function (H3K4me3.ML.PCHIC){
    #Filter out rows where gene is "."
    H3K4me3.ML.PCHIC.genes <- H3K4me3.ML.PCHIC %>% filter (gene !=".")
    #split genes separated by ";" into a single list
    list <- strsplit(H3K4me3.ML.PCHIC.genes$gene, ";")
    # use unlist to put each gene on a new row
    final_list <- as.data.frame(unlist(list))
    #show only unique genes
    unique_list <- unique(final_list)
    #rename columns to "gene"
    colnames(unique_list) <- c("gene")
    return(unique_list)
}

H3K4me3.CD4.background.PCHIC <- extract.H3K4me3.background.PCHIC(H3K4me3.CD4.PCHIC)
H3K4me3.CD8.background.PCHIC <- extract.H3K4me3.background.PCHIC(H3K4me3.CD8.PCHIC)
H3K4me3.CD14.background.PCHIC <- extract.H3K4me3.background.PCHIC(H3K4me3.CD14.PCHIC)
#Make a concatenated background set using rbind
H3K4me3.background.PCHIC <- rbind (H3K4me3.CD4.background.PCHIC, 
                                H3K4me3.CD8.background.PCHIC, 
                                H3K4me3.CD14.background.PCHIC)

#Generate concatenated lists of genes from both methods
H3K4me3.CD4.genes <- rbind (H3K4me3.CD4.Top200.proximity, 
                         H3K4me3.CD4.Top200.PCHIC)
H3K4me3.CD8.genes <- rbind (H3K4me3.CD8.Top200.proximity, 
                         H3K4me3.CD8.Top200.PCHIC)
H3K4me3.CD14.genes <- rbind (H3K4me3.CD14.Top200.proximity, 
                          H3K4me3.CD14.Top200.PCHIC)
#Convert these lists to vectors
H3K4me3.CD4.genes <- as.vector(H3K4me3.CD4.genes$gene)
H3K4me3.CD8.genes <- as.vector(H3K4me3.CD8.genes$gene)
H3K4me3.CD14.genes <- as.vector(H3K4me3.CD14.genes$gene)
length(H3K4me3.CD4.genes) #406
length(H3K4me3.CD8.genes)#372
length(H3K4me3.CD14.genes)#280

#Generate concatenated background for all 3 cell types
H3K4me3.background <- rbind(H3K4me3.background.proximity,
                         H3K4me3.background.PCHIC)
H3K4me3.background <- as.vector(H3K4me3.background$gene)
class(H3K4me3.background)

##Run pathway enrichment using XGR


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

H3K4me3.CD4.GOBP <- XGR_function(H3K4me3.CD4.genes, H3K4me3.background, ontology = "GOBP")
H3K4me3.CD4.GOBP.table <- xEnrichViewer(H3K4me3.CD4.GOBP, details = TRUE, 
                                     top_num = 50)
write.table(H3K4me3.CD4.GOBP.table, "./GOBP/H3K4me3.CD4.GOBP.txt", sep = "\t", quote = FALSE)

H3K4me3.CD8.GOBP <- XGR_function(H3K4me3.CD8.genes, H3K4me3.background, ontology = "GOBP")
H3K4me3.CD8.GOBP.table <- xEnrichViewer(H3K4me3.CD8.GOBP, details = TRUE, 
                                     top_num = 50)
write.table(H3K4me3.CD8.GOBP.table, "./GOBP/H3K4me3.CD8.GOBP.txt", sep = "\t", quote = FALSE)

H3K4me3.CD14.GOBP <- XGR_function(H3K4me3.CD14.genes, H3K4me3.background, ontology = "GOBP")
H3K4me3.CD14.GOBP.table <- xEnrichViewer(H3K4me3.CD14.GOBP, details = TRUE, 
                                      top_num = 50)
write.table(H3K4me3.CD14.GOBP.table, "./GOBP/H3K4me3.CD14.GOBP.txt", sep = "\t", quote = FALSE)

#make a list of things to compare
H3K4me3.GOBP <- list(H3K4me3.CD4.GOBP, H3K4me3.CD8.GOBP, H3K4me3.CD14.GOBP)
names(H3K4me3.GOBP) <- c("CD4 T cells", "CD8 T cells", "CD14 monocytes")

#Set the order
H3K4me3.GOBP <- H3K4me3.GOBP[c("CD14 monocytes", "CD8 T cells", "CD4 T cells")]

H3K4me3.GOBP.compare <- xEnrichCompare(H3K4me3.GOBP, FDR.cutoff = 0.01, facet = TRUE)
H3K4me3.GOBP.compare <- 
    H3K4me3.GOBP.compare + theme(axis.text.y=element_text(size=14,color="black")) +
    scale_fill_manual(values = c(Okabe_Ito[3], Okabe_Ito[5], Okabe_Ito[6])) 
H3K4me3.GOBP.compare

dir.create("figures", showWarnings = FALSE)
ggsave("H3K4me3.GOBP.pdf", H3K4me3.GOBP.compare, width = 15, height = 15, useDingbats = FALSE, path = "./figures/")
ggsave("H3K4me3.GOBP.png", H3K4me3.GOBP.compare, width = 15, height = 15, device = "png", path = "./figures/")


#Now do for Reactome

H3K4me3.CD4.Reactome <- XGR_function(H3K4me3.CD4.genes, H3K4me3.background, ontology = "MsigdbC2REACTOME")

H3K4me3.CD4.Reactome.table <- xEnrichViewer(H3K4me3.CD4.Reactome, details = TRUE, 
                                         top_num = 50)

write.table(H3K4me3.CD4.Reactome.table, "./Reactome/H3K4me3.CD4.Reactome.txt", sep = "\t", quote = FALSE)

H3K4me3.CD8.Reactome <- XGR_function(H3K4me3.CD8.genes, H3K4me3.background, ontology = "MsigdbC2REACTOME")
H3K4me3.CD8.Reactome.table <- xEnrichViewer(H3K4me3.CD8.Reactome, details = TRUE, 
                                         top_num = 50)
write.table(H3K4me3.CD8.Reactome.table, "./Reactome/H3K4me3.CD8.Reactome.txt", sep = "\t", quote = FALSE)

H3K4me3.CD14.Reactome <- XGR_function(H3K4me3.CD14.genes, H3K4me3.background, ontology = "MsigdbC2REACTOME")
H3K4me3.CD14.Reactome.table <- xEnrichViewer(H3K4me3.CD14.Reactome, details = TRUE, 
                                          top_num = 50)
write.table(H3K4me3.CD14.Reactome.table, "./Reactome/H3K4me3.CD14.Reactome.txt", sep = "\t", quote = FALSE)

H3K4me3.Reactome <- list(H3K4me3.CD4.Reactome, H3K4me3.CD8.Reactome, H3K4me3.CD14.Reactome)
names(H3K4me3.Reactome) <- c("CD4+ T cells", "CD8+ T cells", "CD14+ monocytes")
#Set the order
H3K4me3.Reactome <- H3K4me3.Reactome[c("CD14+ monocytes", "CD8+ T cells", "CD4+ T cells")]

#Perform the comparison
H3K4me3.Reactome.compare <- xEnrichCompare(H3K4me3.Reactome, FDR.cutoff = 0.05, 
                                           signature = FALSE,
                                           bar.label.size = 4)

H3K4me3.Reactome.compare <- 
    H3K4me3.Reactome.compare + theme(axis.text.y=element_text(size=18,color="black")) +
    scale_fill_manual(values = c(Okabe_Ito[3], Okabe_Ito[5], Okabe_Ito[6])) 
H3K4me3.Reactome.compare

dir.create("figures", showWarnings = FALSE)
ggsave("H3K4me3.Reactome.pdf", H3K4me3.Reactome.compare, width = 15, height = 10, useDingbats = FALSE, path = "./figures/")
ggsave("H3K4me3.Reactome.png", H3K4me3.Reactome.compare, width = 15, height = 10, device = "png", path = "./figures/")



#Find examples from output 
# e.g. CCR7 is one of two top significant H3K4me3 Differential peaks in CD4
#Is CCR7 in the list of PCHIC genes or proximity genes?
colnames(H3K4me3.CD14.Top200.PCHIC)
H3K4me3.CD14.Top200.PCHIC %>% filter (gene == "CCR7")
#No
H3K4me3.CD14.Top200.proximity %>% filter (gene == "CCR7")
#Yes


#Find location of H3K4me3 peak in proximal set & generate interaction track file of H3K4me3 peaks
CCR7_loc <- H3K4me3.CD14 %>% filter (grepl("CCR7", ProxGene)) %>%
    select(chr, start, end, padj)
CCR7_loc
write.table(CCR7_loc, "./Examples/H3K4me3_CD14_CCR7.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

#Are any Reactome members_Overlap genes in the PCHIC top 200?
Reactome_members_Overlap.CD14 <- str_split(H3K4me3.CD14.Reactome.table$members_Overlap, ",")
Reactome_members_Overlap.CD14 <- unlist(Reactome_members_Overlap.CD14)
Reactome_members_Overlap.CD14 <- as.data.frame(Reactome_members_Overlap.CD14)
colnames(Reactome_members_Overlap.CD14) <- c("gene")
Reactome_members_Overlap.CD14
Reactome_in_PCHIC.CD14 <- Reactome_members_Overlap.CD14 %>% inner_join (H3K4me3.CD14.Top200.PCHIC) %>%
    left_join (H3K4me3.CD14.Top200, by = c("gene" = "ProxGene")) %>%
    select (gene, pvalue, padj) 
Reactome_in_PCHIC.CD14
#Only 2 genes, CCR5 and IFIT1 , both have NS padj


#Find Reactome members overlap that intersect with proximal top 200

Reactome_in_proximal.CD14 <- Reactome_members_Overlap.CD14 %>% 
    inner_join (H3K4me3.CD14.Top200, by = c("gene" = "ProxGene")) %>%
    dplyr::select (gene, pvalue, padj) 
Reactome_in_proximal.CD14
#BTLA, CD247, GRAP2, ANAPC1, AGTR1, CUBN, GNA11 but none have sig padj
#BTLA is the best one padj = 0.0647

H3K4me3.CD14 %>% filter (grepl('BTLA', ProxGene))




###Are any genes from enriched pathways near GWAS loci?### Not using

#Create a function called get_gene_locations
# get locations of genes from members_Overlap using BioMart
# create a df in bed format of the output
library(biomaRt)
mart <- useEnsembl("genes", "hsapiens_gene_ensembl", GRCh = "37")

head(H3K4me3.CD8.GOBP.table[,13])
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
gene_members_list <- list(H3K4me3.CD14.GOBP.table, H3K4me3.CD4.GOBP.table, H3K4me3.CD8.GOBP.table,
                          H3K4me3.CD14.Reactome.table, H3K4me3.CD4.Reactome.table, H3K4me3.CD8.Reactome.table)
gene_members_list_bed <- lapply(gene_members_list, function (x) get_gene_locations(x))
names(gene_members_list_bed) <- c("H3K4me3.CD14.GOBP", "H3K4me3.CD4.GOBP", "H3K4me3.CD8.GOBP",
                                  "H3K4me3.CD14.Reactome", "H3K4me3.CD4.Reactome", "H3K4me3.CD8.Reactome")
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

###Are any regions near GWAS loci?### Another method
#Used bedtools to look for any H3K4me3 peaks within 500kb of a lead SNP
#Results are in ./peak_vs_SNP

#Read in the data

H3K4me3_ML_CD4.SNP <- read.table("./peak_vs_SNP/H3K4me3_ML_CD4.SNP.txt", header = FALSE, sep = "\t")
colnames (H3K4me3_ML_CD4.SNP) <- c("chr.SNP", "start.SNP", "end.SNP", "SNP", 
                                "chr", "start", "end", "1")
H3K4me3_ML_CD8.SNP <- read.table("./peak_vs_SNP/H3K4me3_ML_CD8.SNP.txt", header = FALSE, sep = "\t")
colnames (H3K4me3_ML_CD8.SNP) <- c("chr.SNP", "start.SNP", "end.SNP", "SNP", 
                                "chr", "start", "end", "1")
H3K4me3_ML_CD14.SNP <- read.table("./peak_vs_SNP/H3K4me3_ML_CD14.SNP.txt", header = FALSE, sep = "\t")
colnames (H3K4me3_ML_CD14.SNP) <- c("chr.SNP", "start.SNP", "end.SNP", "SNP", 
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

H3K4me3.CD4.SNP <- SNPs_peaks(H3K4me3_ML_CD4.SNP, H3K4me3.CD4)
H3K4me3.CD8.SNP <- SNPs_peaks(H3K4me3_ML_CD8.SNP, H3K4me3.CD8)
H3K4me3.CD14.SNP <- SNPs_peaks(H3K4me3_ML_CD14.SNP, H3K4me3.CD14)
dim(H3K4me3.CD4.SNP)
dim(H3K4me3.CD8.SNP)
dim(H3K4me3.CD14.SNP)

#Lots of interesting examples to check out
#CD14 IL7R at promoter
#POFUT1 in gene body
#NUPR1 not great
#SPHK2 v low peak in gene body
#LY9 at promoter but far from cred set


write.table (H3K4me3.CD4.SNP, "./peak_vs_SNP/H3K4me3.CD4.SNP.txt", 
             col.names = TRUE,
             row.names = FALSE,
             quote = FALSE,
             na = "NA",
             sep = "\t")

write.table (H3K4me3.CD8.SNP, "./peak_vs_SNP/H3K4me3.CD8.SNP.txt", 
             col.names = TRUE,
             row.names = FALSE,
             quote = FALSE,
             na = "NA",
             sep = "\t")
write.table (H3K4me3.CD14.SNP, "./peak_vs_SNP/H3K4me3.CD14.SNP.txt", 
             col.names = TRUE,
             row.names = FALSE,
             quote = FALSE,
             na = "NA",
             sep = "\t")

#Are any H3K4me3 differential peaks at DE genes?
#Read in gene expression data

#Import DEseq2 data from RNA-seq

RNA.CD4 <- read.table("C:/Users/ccohen.WHGC/Documents/R/XGR/RNAseq_cohort2-5/Inputs/RNA.CD4.txt",
                      header = TRUE, sep = "\t") 
RNA.CD8 <- read.table("C:/Users/ccohen.WHGC/Documents/R/XGR/RNAseq_cohort2-5/Inputs/RNA.CD8.txt",
                      header = TRUE, sep = "\t") 
RNA.CD14 <- read.table("C:/Users/ccohen.WHGC/Documents/R/XGR/RNAseq_cohort2-5/Inputs/RNA.CD14.txt",
                       header = TRUE, sep = "\t") 
#Intersect Top 200 H3K4me3 peaks with RNAseq data and select RNA >1.5 fold change and padj<0.05
CD4.H3K4me3.RNA.sig <- H3K4me3.CD4.Top200 %>% inner_join(RNA.CD4, by = c("ProxGene" = "name")) %>%
    filter(log2FoldChange.y > 0.585 | log2FoldChange.y < -0.585) %>%
    filter(padj.y<0.05) 
CD8.H3K4me3.RNA.sig <- H3K4me3.CD8.Top200 %>% inner_join(RNA.CD8, by = c("ProxGene" = "name")) %>%
    filter(log2FoldChange.y > 0.585 | log2FoldChange.y < -0.585) %>%
    filter(padj.y<0.05) 
CD14.H3K4me3.RNA.sig <- H3K4me3.CD14.Top200 %>% inner_join(RNA.CD14, by = c("ProxGene" = "name")) %>%
    filter(log2FoldChange.y > 0.585 | log2FoldChange.y < -0.585) %>%
    filter(padj.y<0.05) 
#CD4 and CD8 have no overlap
#CD14 has 12 regions
dim(CD14.H3K4me3.RNA.sig)
head(CD14.H3K4me3.RNA.sig)

write.table (CD14.H3K4me3.RNA.sig, "./Examples/CD14.H3K4me3.RNA.sig.txt", 
             col.names = TRUE,
             row.names = FALSE,
             quote = FALSE,
             na = "NA",
             sep = "\t")
