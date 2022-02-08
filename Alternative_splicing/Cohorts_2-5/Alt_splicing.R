#Look at alt splicing results from cohorts 2-5

library(tidyverse)

# read the files
library(data.table)
Alt_splicing.CD4 <- fread("./DEXseq/DEXSeq.results.CD4_0721.csv")
class (Alt_splicing.CD4$genomicData.end)
Alt_splicing.CD8 <- fread("./DEXseq/DEXSeq.results.CD8_0721.csv")
class (Alt_splicing.CD8$genomicData.end)
Alt_splicing.CD14 <- fread("./DEXseq/DEXSeq.results.CD14_0721.csv")
Alt_splicing.CD19 <- fread("./DEXseq/DEXSeq.results.CD19_0721.csv")

Alt_splicing <- list(Alt_splicing.CD4, Alt_splicing.CD8, Alt_splicing.CD14, Alt_splicing.CD19)
names(Alt_splicing) <- c("CD4", "CD8", "CD14", "CD19")

#select significant exons
select_sig_exons <- function (df) {
    df2 <- df %>% drop_na(padj, log2fold_AS_HV) %>%
        filter(padj < 0.05) %>%
        filter(log2fold_AS_HV < -0.585 | log2fold_AS_HV > 0.585) %>%
        filter (genomicData.strand != "") %>%
        select (1:16)
    return (df2)
}

Alt_splicing.sig <- lapply(Alt_splicing, select_sig_exons)

Alt_splicing.CD4.sig <- as.data.frame(Alt_splicing.sig[1])

#How many significant regions are there?
dim(Alt_splicing.sig[[1]]) #CD4 315
dim(Alt_splicing.sig[[2]]) #CD8 371
dim(Alt_splicing.sig[[3]]) #CD14 1286
dim(Alt_splicing.sig[[4]]) #CD19 0


#Convert EnsemblID to gene name using BiomaRt
library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
#select ensemblIDs from the list of Alt splicing sig results
ensemblIDs <- lapply(Alt_splicing.sig, function (x) dplyr::select(x, 2))

gene_IDs <- list()
for (i in 1:(length(Alt_splicing.sig))) {
    gene_IDs[[i]] <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                  values = ensemblIDs[[i]], mart= mart)
}

#NB get an error message because the list of ensemblIDs in CD19 is length 0

#Rename column in gene_Ids
gene_IDs <- lapply(gene_IDs, function (x) dplyr::rename (x, groupID = ensembl_gene_id))

Alt_splicing.sig.gene <- list()
for (i in 1:(length(Alt_splicing.sig))) {
    Alt_splicing.sig.gene[[i]] <- Alt_splicing.sig[[i]] %>% left_join (gene_IDs[[i]])
}
#save results

dir.create("Results")
for (i in 1:(length(Alt_splicing.sig))) {
    write.table (Alt_splicing.sig.gene[i], 
                 paste0("./Results/Alt_splicing.sig.", names(Alt_splicing.sig[i]), ".txt"),
                 col.names = TRUE,
                 row.names = FALSE,
                 quote = FALSE, 
                 sep = "\t")
}


#Create a bed file so I can use bedtools to check if there are any diff alt splicing events in GWAS regions

make_bed <- function (df) {
    df2 <- df %>% 
        mutate (chr = paste0("chr", genomicData.seqnames)) %>%
        dplyr::select (chr, genomicData.start, genomicData.end, hgnc_symbol) %>%
        mutate_all(na_if,"") %>%
        #remove rows where col 3 == col 2
        filter (genomicData.start != genomicData.end)
    return (df2)
}

Alt_splicing.sig.bed <- lapply(Alt_splicing.sig.gene, make_bed)


#save bed files

dir.create("bed")
for (i in 1:(length(Alt_splicing.sig))) {
    write.table (Alt_splicing.sig.bed[i], 
                 paste0("./bed/Alt_splicing.sig.", names(Alt_splicing.sig[i]), ".bed"),
                 col.names = FALSE,
                 row.names = FALSE,
                 quote = FALSE,
                 na = "NA",
                 sep = "\t")
}


Alt_splicing.sig.CD14 <- as.data.frame(Alt_splicing.sig.gene[[3]])


#use 
# >  bedtools window -a LeadSNPs.sorted.bed -b Alt_splicing.sig.CD14.sort.bed -w 100000 > ./window_100000/Alt_splicing.sig.CD14_vs_leadSNPs.txt
#etc

#saved in folder Alt_splicing_vs_GWAS
#window_100000 and window_50000 for window size

CD4.GWAS.100000 <- read.table ("./Alt_splicing_vs_GWAS/window_100000/Alt_splicing.sig.CD4_vs_leadSNPs.txt", 
            sep = "\t")
CD8.GWAS.100000 <- read.table ("./Alt_splicing_vs_GWAS/window_100000/Alt_splicing.sig.CD8_vs_leadSNPs.txt", 
                                sep = "\t")
CD14.GWAS.100000 <- read.table ("./Alt_splicing_vs_GWAS/window_100000/Alt_splicing.sig.CD14_vs_leadSNPs.txt", 
                                sep = "\t")
CD4.GWAS.500000 <- read.table ("./Alt_splicing_vs_GWAS/window_500000/Alt_splicing.sig.CD4_vs_leadSNPs.txt", 
                               sep = "\t")
CD8.GWAS.500000 <- read.table ("./Alt_splicing_vs_GWAS/window_500000/Alt_splicing.sig.CD8_vs_leadSNPs.txt", 
                               sep = "\t")
CD14.GWAS.500000 <- read.table ("./Alt_splicing_vs_GWAS/window_500000/Alt_splicing.sig.CD14_vs_leadSNPs.txt", 
                                sep = "\t")

Alt_splicing.GWAS.500000 <- list (CD4.GWAS.500000, CD8.GWAS.500000, CD14.GWAS.500000)

#How many significant alt_splicing events are within 500kb of a lead SNP?
dim(CD4.GWAS.500000) #CD4 21
dim(CD8.GWAS.500000) #CD8 22
dim(CD14.GWAS.500000) #CD14 87

column_names <- c("chr.SNP", "start.SNP", "end.SNP", "SNP", "chr", "genomicData.start", "genomicData.end", "gene")


Alt_splicing.GWAS.500000 <- lapply(Alt_splicing.GWAS.500000, setNames, column_names)

#Join regions from GWAS with full results
head(Alt_splicing.GWAS.500000[1])
head(Alt_splicing.sig[1])

Alt_splicing.CD4$genomicData.start <- as.integer(Alt_splicing.CD4$genomicData.start)
Alt_splicing.CD4$genomicData.end <- as.integer(Alt_splicing.CD4$genomicData.end)
Alt_splicing.CD8$genomicData.start <- as.integer(Alt_splicing.CD8$genomicData.start)
Alt_splicing.CD8$genomicData.end <- as.integer(Alt_splicing.CD8$genomicData.end)
Alt_splicing.CD14$genomicData.start <- as.integer(Alt_splicing.CD14$genomicData.start)
Alt_splicing.CD14$genomicData.end <- as.integer(Alt_splicing.CD14$genomicData.end)

Alt_splicing.GWAS.500000.CD4 <- Alt_splicing.GWAS.500000[[1]] %>% left_join (Alt_splicing.CD4) %>%
    select (1:22)
Alt_splicing.GWAS.500000.CD8 <- Alt_splicing.GWAS.500000[[2]] %>% left_join (Alt_splicing.CD8) %>%
    select (1:22)
Alt_splicing.GWAS.500000.CD14 <- Alt_splicing.GWAS.500000[[3]] %>% left_join (Alt_splicing.CD14) %>%
    select (1:22)

head(Alt_splicing.GWAS.500000.CD4)
write.table(Alt_splicing.GWAS.500000.CD4, 
            "./Alt_splicing_vs_GWAS/window_500000/Alt_splicing.GWAS.5000000.CD4.txt", 
            quote = FALSE, 
            row.names = FALSE, 
            sep = "\t")
write.table(Alt_splicing.GWAS.500000.CD8, 
            "./Alt_splicing_vs_GWAS/window_500000/Alt_splicing.GWAS.5000000.CD8.txt", 
            quote = FALSE, 
            row.names = FALSE, 
            sep = "\t")
write.table(Alt_splicing.GWAS.500000.CD14, 
            "./Alt_splicing_vs_GWAS/window_500000/Alt_splicing.GWAS.5000000.CD14.txt", 
            quote = FALSE, 
            row.names = FALSE, 
            sep = "\t")

