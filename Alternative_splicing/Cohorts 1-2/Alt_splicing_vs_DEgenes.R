#First way round
#In Ping's existing alt splicing significant results, are there any DE genes?

#Read in alt splicing results
Alt_splicing_sig.CD14 <- read.table("Alt_splicing_sig.CD14.txt", header = TRUE, sep = "\t")
#rename column
library(dplyr)
Alt_splicing_sig.CD14 <- Alt_splicing_sig.CD14 %>% rename(name = hgnc_symbol)

#Read in DE genes
DEgenes.CD14 <- read.table("DEgenes.CD14.txt", header = TRUE, sep = "\t")

#See if there are any genes in both tables
merge.CD14 <- DEgenes.CD14 %>% inner_join(Alt_splicing_sig.CD14, by = c("name"))

#make bed file of significant regions
merge.CD14.bed <- merge.CD14 %>% select(genomicData.seqnames, genomicData.start, genomicData.end, name)
merge.CD14.bed$genomicData.seqnames <- paste("chr", merge.CD14.bed$genomicData.seqnames, sep = "")
colnames(merge.CD14.bed) <- c("chr", "start", "end", "name")
write.table(merge.CD14.bed, "merge.CD14.bed", col.names = FALSE, row.names=FALSE,sep="\t", quote = FALSE)

#Second way round

#Use bash script to select all DE genes in Ping's complete data set which is
# /well/jknight/ping/AS.patinets.splicing/DEXSeq.results.CD14_cohort1.2.csv
#(I have made a softlink to this file in my folder /well/jknight/Carla/Alt_splicing)

#Read in these results
Alt_splicing_DEgenes.CD14 <- read.table("Alt_splicing_DEgenes.CD14.txt", header = FALSE, sep = "\t")
#rename the columns
colnames(Alt_splicing_DEgenes.CD14) <- c("name","groupID", 
                                            "featureID",
                                            "exonBaseMean",
                                            "dispersion",
                                            "stat",
                                            "pvalue",
                                            "padj",
                                            "controls",
                                            "cases",
                                            "log2fold_cases_controls",
                                            "genomicData.seqnames",
                                            "genomicData.start",
                                            "genomicData.end",
                                            "genomicData.width",
                                            "genomicData.strand",
                                            "countData.22",
                                            "countData.23",
                                            "countData.1",
                                            "countData.2",
                                            "countData.3",
                                            "countData.4",
                                            "countData.6",
                                            "countData.7",
                                            "countData.8",
                                            "countData.9",
                                            "countData.10",
                                            "countData.12",
                                            "countData.13",
                                            "countData.14",
                                            "countData.15",
                                            "countData.16",
                                            "countData.17",
                                            "countData.18",
                                            "countData.19",
                                            "countData.20",
                                            "countData.21",
                                            "countData.24",
                                            "countData.25", 
                                            "countData.26",
                                            "countData.27",
                                            "countData.28",
                                            "countData.29",
                                            "countData.30",
                                            "countData.31",
                                            "countData.32",
                                            "countData.33",
                                            "countData.34",
                                            "countData.35",
                                            "countData.36",
                                            "countData.37",
                                            "countData.38",
                                            "countData.39",
                                            "countData.40", 
                                            "countData.41",
                                            "countData.42",
                                            "countData.43",
                                            "countData.44",
                                            "countData.45",
                                            "countData.46",
                                            "countData.47",
                                            "countData.48",
                                            "countData.49",
                                            "countData.50",
                                            "countData.51",
                                            "countData.52",
                                            "countData.53",
                                            "countData.54")

#select columns 1-69
Alt_splicing_DEgenes.CD14 <- Alt_splicing_DEgenes.CD14[,1:68]
#Convert EnsemblID to gene name using BiomaRt
library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
ensemblIDs <-  Alt_splicing_DEgenes.CD14$groupID
gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                  values = ensemblIDs, mart= mart)
library(dplyr)
gene_IDs <- gene_IDs %>% dplyr::rename (groupID = ensembl_gene_id)

Alt_splicing_DEgenes.CD14 <- gene_IDs %>% full_join(Alt_splicing_DEgenes.CD14)
Alt_splicing_DEgenes.CD14.p001.FC1.5 <- Alt_splicing_DEgenes.CD14 %>% 
    filter(pvalue < 0.01) %>% 
    filter(log2fold_cases_controls > 0.585 | log2fold_cases_controls < -0.585)
#save filtered results
write.table(Alt_splicing_DEgenes.CD14.p001.FC1.5, "Alt_splicing_DEgenes.CD14.p001.FC1.5", row.names = FALSE, sep = "\t", quote = FALSE)


#Next approach is to identify all exons near GWAS regions (done in bedtools)
Alt_splicing_near_SNPs.CD14 <- read.table("Alt_splicing_near_SNPs.csv", header = TRUE, sep = ",")
#Need to change EnsemblID to gene name
library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
ensemblIDs_near_SNPs <-  Alt_splicing_near_SNPs.CD14$exonID
gene_IDs_SNPs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                  values = ensemblIDs_near_SNPs, mart= mart)
library(dplyr)
gene_IDs_SNPs <- gene_IDs_SNPs %>% dplyr::rename (exonID = ensembl_gene_id)
Alt_splicing_near_SNPs.CD14 <- gene_IDs_SNPs %>% full_join(Alt_splicing_near_SNPs.CD14)
Alt_splicing_near_SNPs.CD14.p001.FC1.5 <- Alt_splicing_near_SNPs.CD14 %>% 
    filter(pvalue < 0.01) %>% 
    filter(log2fold_cases_controls > 0.585 | log2fold_cases_controls < -0.585)

#Filter log2FC only
Alt_splicing_near_SNPs.CD14.FC1.5 <- Alt_splicing_near_SNPs.CD14 %>% 
        filter(log2fold_cases_controls > 0.585 | log2fold_cases_controls < -0.585)

#save filtered results
write.table(Alt_splicing_near_SNPs.CD14.p001.FC1.5, "Alt_splicing_near_SNPs.CD14.p001.FC1.5", row.names = FALSE, sep = "\t", quote = FALSE)

#Make bed file of CCRL2 results
#make bed file of significant regions
CCRL2_bed  <- Alt_splicing_near_SNPs.CD14.p001.FC1.5 %>% 
    select(genomicData.seqnames, genomicData.start, genomicData.end, hgnc_symbol) %>%
    filter (hgnc_symbol == "CCRL2")
CCRL2_bed$genomicData.seqnames <- paste("chr", CCRL2_bed$genomicData.seqnames, sep = "")
colnames(CCRL2_bed) <- c("chr", "start", "end", "name")
write.table(CCRL2_bed, "CCRL2.CD14.bed", col.names = FALSE, row.names=FALSE,sep="\t", quote = FALSE)

