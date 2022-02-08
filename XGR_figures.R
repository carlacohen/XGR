#Aim to draw graph for XGR output
setwd("~/R/XGR/")
library(tidyverse)
library(ggrepel)
library(cowplot)

#Define colours
Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000") 
#CD4 009E73 Okabe_Ito[3] blueish green
#CD8 0072B2 Okabe_Ito[5] blue
#CD14 D55E00 Okabe_Ito[6] vermillion

#Read in files
ATAC.CD4.Reactome <- read.table("./ATAC_cohort2-5/Reactome/ATAC.CD4.Reactome_2.txt",
                                header = TRUE, sep = "\t")
ATAC.CD8.Reactome <- read.table("./ATAC_cohort2-5/Reactome/ATAC.CD8.Reactome_2.txt",
                                header = TRUE, sep = "\t")
ATAC.CD14.Reactome <- read.table("./ATAC_cohort2-5/Reactome/ATAC.CD14.Reactome_2.txt",
                                 header = TRUE, sep = "\t")
H3K4me3.CD4.Reactome <- read.table("./H3K4me3_cohort2-5/Reactome/H3K4me3.CD4.Reactome.txt",
                                header = TRUE, sep = "\t")
H3K4me3.CD8.Reactome <- read.table("./H3K4me3_cohort2-5/Reactome/H3K4me3.CD8.Reactome.txt",
                                header = TRUE, sep = "\t")
H3K4me3.CD14.Reactome <- read.table("./H3K4me3_cohort2-5/Reactome/H3K4me3.CD14.Reactome.txt",
                                 header = TRUE, sep = "\t")
H3K27ac.CD4.Reactome <- read.table("./H3K27ac_cohort2-5/Reactome/H3K27ac.CD4.Reactome.txt",
                                   header = TRUE, sep = "\t")
H3K27ac.CD8.Reactome <- read.table("./H3K27ac_cohort2-5/Reactome/H3K27ac.CD8.Reactome.txt",
                                   header = TRUE, sep = "\t")
H3K27ac.CD14.Reactome <- read.table("./H3K27ac_cohort2-5/Reactome/H3K27ac.CD14.Reactome.txt",
                                    header = TRUE, sep = "\t")
eRNA.CD4.Reactome <- read.table("./eRNA_cohort2-5/Reactome/eRNA.CD4.Reactome.txt",
                                   header = TRUE, sep = "\t")
eRNA.CD8.Reactome <- read.table("./eRNA_cohort2-5/Reactome/eRNA.CD8.Reactome.txt",
                                   header = TRUE, sep = "\t")
eRNA.CD14.Reactome <- read.table("./eRNA_cohort2-5/Reactome/eRNA.CD14.Reactome.txt",
                                    header = TRUE, sep = "\t")

#Make a combined df in order to set the levels of the named categories

make_df1 <- function (df, celltype) {
    df1 <- df %>%
        filter (adjp < 0.01) %>%
        arrange (fc) %>%
        mutate ("{celltype}.log2fc" := log2(fc), .keep = "unused") %>% 
        mutate ("{celltype}.padj" := log10(adjp) * -1, .keep = "unused") %>%
        mutate (cell_type := {{celltype}}, .keep = "unused") %>%
        rename ("{celltype}.nOverlap" := nOverlap) %>%
        rename ("{celltype}.zscore" := zscore) %>%
        select (!c(nAnno, pvalue, or, CIl, CIu, distance, namespace, members_Overlap, members_Anno, cell_type))
    return(df1)
}

ATAC.CD4.Reactome.select <- make_df1(ATAC.CD4.Reactome, "CD4")
ATAC.CD8.Reactome.select <- make_df1(ATAC.CD8.Reactome, "CD8")
ATAC.CD14.Reactome.select <- make_df1(ATAC.CD14.Reactome, "CD14")
H3K4me3.CD4.Reactome.select <- make_df1(H3K4me3.CD4.Reactome, "CD4")
H3K4me3.CD8.Reactome.select <- make_df1(H3K4me3.CD8.Reactome, "CD8")
H3K4me3.CD14.Reactome.select <- make_df1(H3K4me3.CD14.Reactome, "CD14")
H3K27ac.CD4.Reactome.select <- make_df1(H3K27ac.CD4.Reactome, "CD4")
H3K27ac.CD8.Reactome.select <- make_df1(H3K27ac.CD8.Reactome, "CD8")
H3K27ac.CD14.Reactome.select <- make_df1(H3K27ac.CD14.Reactome, "CD14")
eRNA.CD4.Reactome.select <- make_df1(eRNA.CD4.Reactome, "CD4")
eRNA.CD8.Reactome.select <- make_df1(eRNA.CD8.Reactome, "CD8")
eRNA.CD14.Reactome.select <- make_df1(eRNA.CD14.Reactome, "CD14")

join_fn <- function (df.CD4, df.CD8, df.CD14) {
    df1 <- df.CD4 %>% 
        full_join (df.CD8) %>% 
        full_join (df.CD14) %>%
        #add a column called "group" to show how many cell types that pathway is present in
        mutate (
            group = case_when(
                CD4.log2fc > 0 & CD8.log2fc > 0 & CD14.log2fc > 0 ~ "1_CD4_CD8_CD14",
                CD8.log2fc > 0 & CD14.log2fc > 0 ~ "2_CD8_CD14",
                CD4.log2fc > 0 & CD14.log2fc > 0 ~ "3_CD4_CD14",
                CD4.log2fc > 0 & CD8.log2fc > 0 ~ "4_CD4_CD8",
                CD14.log2fc > 0 ~ "5_CD14", 
                CD8.log2fc > 0 ~ "6_CD8",
                CD4.log2fc > 0 ~ "7_CD4"
            )
        )  %>% 
        arrange (desc(CD4.log2fc), desc(CD8.log2fc), desc(CD14.log2fc)) %>%
        arrange (group)
    #Make "name" an ordered factor
    df1$name <- factor (df1$name, levels = df1$name)
    return (df1)
}

ATAC.Reactome <- join_fn(ATAC.CD4.Reactome.select, ATAC.CD8.Reactome.select, ATAC.CD14.Reactome.select)
H3K4me3.Reactome <- join_fn(H3K4me3.CD4.Reactome.select, H3K4me3.CD8.Reactome.select, H3K4me3.CD14.Reactome.select)
H3K27ac.Reactome <- join_fn(H3K27ac.CD4.Reactome.select, H3K27ac.CD8.Reactome.select, H3K27ac.CD14.Reactome.select)
eRNA.Reactome <- join_fn(eRNA.CD4.Reactome.select, eRNA.CD8.Reactome.select, eRNA.CD14.Reactome.select)

#Make a second df for plotting, organise data by log2fc and FDR
make_df2 <- function (df, celltype){
    df1 <- df %>% 
        filter (adjp < 0.01) %>%
        mutate (log2fc = log2(fc)) %>% 
        mutate (minuslog10padj = log10(adjp) * -1) %>%
        mutate (cell_type := {{celltype}}) %>%
        select (name, log2fc, minuslog10padj, nOverlap, cell_type, zscore) %>%
        arrange (log2fc)
    
    
}

ATAC.CD4.Reactome.select2 <- make_df2(ATAC.CD4.Reactome, "CD4")
ATAC.CD8.Reactome.select2 <- make_df2(ATAC.CD8.Reactome, "CD8")
ATAC.CD14.Reactome.select2 <- make_df2(ATAC.CD14.Reactome, "CD14")
H3K4me3.CD4.Reactome.select2 <- make_df2(H3K4me3.CD4.Reactome, "CD4")
H3K4me3.CD8.Reactome.select2 <- make_df2(H3K4me3.CD8.Reactome, "CD8")
H3K4me3.CD14.Reactome.select2 <- make_df2(H3K4me3.CD14.Reactome, "CD14")
H3K27ac.CD4.Reactome.select2 <- make_df2(H3K27ac.CD4.Reactome, "CD4")
H3K27ac.CD8.Reactome.select2 <- make_df2(H3K27ac.CD8.Reactome, "CD8")
H3K27ac.CD14.Reactome.select2 <- make_df2(H3K27ac.CD14.Reactome, "CD14")
eRNA.CD4.Reactome.select2 <- make_df2(eRNA.CD4.Reactome, "CD4")
eRNA.CD8.Reactome.select2 <- make_df2(eRNA.CD8.Reactome, "CD8")
eRNA.CD14.Reactome.select2 <- make_df2(eRNA.CD14.Reactome, "CD14")

join_factor_fn <- function (df.CD4, df.CD8, df.CD14, df.join) {
    df1 <- df.CD4 %>% 
        full_join(df.CD8) %>%
        full_join(df.CD14)
    #Set the levels of cell_type and name to get the correct display order
    df1$cell_type <- factor(df1$cell_type, levels = c("CD4", "CD8", "CD14"))
    df1$name <- factor(df1$name, levels = df.join$name) 
    #Reverse the factor levels of name
    df1$name <- fct_rev(df1$name)
    return(df1)
}

ATAC.Reactome2 <- join_factor_fn(ATAC.CD4.Reactome.select2, ATAC.CD8.Reactome.select2, 
                                 ATAC.CD14.Reactome.select2, ATAC.Reactome)
H3K4me3.Reactome2 <- join_factor_fn(H3K4me3.CD4.Reactome.select2, H3K4me3.CD8.Reactome.select2, 
                                 H3K4me3.CD14.Reactome.select2, H3K4me3.Reactome)
H3K27ac.Reactome2 <- join_factor_fn(H3K27ac.CD4.Reactome.select2, H3K27ac.CD8.Reactome.select2, 
                                 H3K27ac.CD14.Reactome.select2, H3K27ac.Reactome)
eRNA.Reactome2 <- join_factor_fn(eRNA.CD4.Reactome.select2, eRNA.CD8.Reactome.select2, 
                                 eRNA.CD14.Reactome.select2, eRNA.Reactome)

#Make a function to draw plots of name vs log2fc with size as -log10padj
plot_fn <- function (df, scale){
    plot <- ggplot(df, aes(x = name, y = log2fc, size = minuslog10padj))+
        geom_point(aes(colour = cell_type)) +
        #scale_size(name = "-log(10)padj", range = c(1, 10))+
        scale_size_continuous(range = c(1,15), breaks = c(5, 10, 15)) + 
        scale_y_continuous(name = "log2(FC)", limits = scale) +
        scale_colour_manual(values = c(CD4 = Okabe_Ito[3], CD8 = Okabe_Ito[5], CD14 = Okabe_Ito[6])) +
        coord_flip() +
        theme_classic() + 
        facet_wrap(~cell_type, drop = FALSE) + 
        theme(#axis.text.y=element_text(size=18,color="black"), 
              #axis.text.x = element_text(size = 18, colour = "black"), 
              #legend.title = element_text(size = 18, colour = "black"), 
              #axis.title = element_text(size = 18, colour = "black"),
              #legend.text = element_text(size = 18, colour = "black"),
              axis.title.y = element_blank(),
              strip.text = element_blank())
    return (plot)
}

#Run the function on each dataset
ATAC.Reactome.plot.padj <- plot_fn(ATAC.Reactome2, scale = c(0, 6))
H3K4me3.Reactome.plot.padj <- plot_fn(H3K4me3.Reactome2, scale = c(0, 5))
H3K27ac.Reactome.plot.padj <- plot_fn(H3K27ac.Reactome2, scale = c(0, 3))
eRNA.Reactome.plot.padj <- plot_fn(eRNA.Reactome2, scale = c(0, 3))

#Combine the 4 plots into a figure using cowplot
figure <- plot_grid(ATAC.Reactome.plot.padj, H3K4me3.Reactome.plot.padj, 
                    H3K27ac.Reactome.plot.padj, eRNA.Reactome.plot.padj,
                    align = "v", 
                    labels = c("ATAC", "H3K4me3", "H3K27Ac", "eRNA"))
figure


ggsave("figure.Reactome.padj.pdf", figure, width = 15, height = 6, useDingbats = FALSE, path = "./Paper_figures/")
ggsave("figure.Reactome.padj.png", figure, width = 15, height = 6, device = "png", path = "./Paper_figures/")


#Save individual plots
dir.create("Paper_figures")
ggsave("ATAC.Reactome.padj.pdf", ATAC.Reactome.plot.padj, width = 15, height = 6, useDingbats = FALSE, path = "./Paper_figures/")
ggsave("ATAC.Reactome.padj.png", ATAC.Reactome.plot.padj, width = 15, height = 6, device = "png", path = "./Paper_figures/")
ggsave("H3K4me3.Reactome.padj.pdf", H3K4me3.Reactome.plot.padj, width = 15, height = 6, useDingbats = FALSE, path = "./Paper_figures/")
ggsave("H3K4me3.Reactome.padj.png", H3K4me3.Reactome.plot.padj, width = 15, height = 6, device = "png", path = "./Paper_figures/")
ggsave("H3K27ac.Reactome.padj.pdf", H3K27ac.Reactome.plot.padj, width = 15, height = 6, useDingbats = FALSE, path = "./Paper_figures/")
ggsave("H3K27ac.Reactome.padj.png", H3K27ac.Reactome.plot.padj, width = 15, height = 6, device = "png", path = "./Paper_figures/")
ggsave("eRNA.Reactome.padj.pdf", eRNA.Reactome.plot.padj, width = 15, height = 6, useDingbats = FALSE, path = "./Paper_figures/")
ggsave("eRNA.Reactome.padj.png", eRNA.Reactome.plot.padj, width = 15, height = 6, device = "png", path = "./Paper_figures/")


##Attempt to make the size of the dots consistent between plots 
#by putting all data into one df and using facet wrap

#First add a column to the Reactome2 table with the method name

ATAC.Reactome3 <- ATAC.Reactome2 %>% mutate (method = "ATAC") %>% arrange (minuslog10padj)
H3K4me3.Reactome3 <- H3K4me3.Reactome2 %>% mutate (method = "H3K4me3") %>% arrange (minuslog10padj)
H3K27ac.Reactome3 <- H3K27ac.Reactome2 %>% mutate (method = "H3K27ac")%>% arrange (minuslog10padj)
eRNA.Reactome3 <- eRNA.Reactome2 %>% mutate (method = "eRNA")%>% arrange (minuslog10padj)

Reactome3.combined <- ATAC.Reactome3 %>% full_join (H3K4me3.Reactome3) %>%
    full_join(H3K27ac.Reactome3) %>%
    full_join(eRNA.Reactome3)

#make the method list a factor with levels to set the order
Reactome3.combined$method <- factor (Reactome3.combined$method, levels = c("ATAC", "H3K4me3", "H3K27ac", "eRNA"))

#plot log2FC, size is minuslog10padj
combined_figure <- ggplot(Reactome3.combined, aes(x = name, y = log2fc, size = minuslog10padj))+
    geom_point(aes(colour = cell_type)) +
    scale_size(name = "-log(10)padj", range = c(1, 10))+
    scale_y_continuous(name = "log2(FC)", limits = c(0,6)) +
    scale_colour_manual(values = c(CD4 = Okabe_Ito[3], CD8 = Okabe_Ito[5], CD14 = Okabe_Ito[6])) +
    coord_flip() +
    theme_classic() +
    facet_grid(. ~ method + cell_type, drop = FALSE) +
    theme(axis.text.y=element_text(size=14,color="black"), 
        axis.text.x = element_text(size = 18, colour = "black"), 
        legend.title = element_text(size = 18, colour = "black"), 
        axis.title = element_text(size = 18, colour = "black"),
        legend.text = element_text(size = 18, colour = "black"),
        axis.title.y = element_blank())

combined_figure
combined_figure + theme (panel.spacing.x = unit(0.5, "lines"))
    
ggsave("combined_figure.Reactome.padj.pdf", combined_figure, width = 18, height = 7, useDingbats = FALSE, path = "./Paper_figures/")
ggsave("combined_figure.Reactome.padj.png", combined_figure, width = 18, height = 7, device = "png", path = "./Paper_figures/")


#Draw a plot of FDR vs z-score with size by number of genes

Reactome.plot.zscore <- ggplot(Reactome3.combined, aes(x = zscore, y = minuslog10padj, size = nOverlap)) + 
    geom_point(aes(colour = cell_type)) +
    scale_size(name = "-log(10)padj", range = c(1, 10))+
    scale_colour_manual(values = c(CD4 = Okabe_Ito[3], CD8 = Okabe_Ito[5], CD14 = Okabe_Ito[6])) +
    scale_y_continuous(name = "-log10(padj)") +
    scale_x_continuous(name = "Z-score") +
    theme_classic() +
    facet_wrap(~ method, scales = "free", nrow = 1) +
    geom_text_repel(
        aes(label = name),
        color = "black",
        size = 10/.pt, # font size 9 pt
        point.padding = 0.1, 
        box.padding = 0.6,
        max.overlaps = 20
    )+
    theme(axis.text.y=element_text(size=24,color="black"), 
          axis.text.x = element_text(size = 24, colour = "black"), 
          legend.title = element_text(size = 24, colour = "black"), 
          axis.title = element_text(size = 24, colour = "black"),
          legend.text = element_text(size = 24, colour = "black"))
Reactome.plot.zscore

ggsave("combined.Reactome.zscore.pdf", Reactome.plot.zscore, width = 20, height = 5, useDingbats = FALSE, path = "./Paper_figures/")
ggsave("combined.Reactome.zscore.png", Reactome.plot.zscore, width = 12, height = 10, device = "png", path = "./Paper_figures/")


#Draw plots of name vs -log10padj with size as nOverlap

combined_figure.padj <- ggplot(Reactome3.combined, aes(x = name, y = minuslog10padj, size = nOverlap))+
    geom_point(aes(colour = cell_type)) +
    scale_size(name = "nOverlap", range = c(1, 10))+
    scale_y_continuous(name = "-log(10)padj", limits = c(0,20), breaks = c(0, 10, 20), labels = c("", "10", "20")) +
    scale_colour_manual(values = c(CD4 = Okabe_Ito[3], CD8 = Okabe_Ito[5], CD14 = Okabe_Ito[6])) +
    coord_flip() +
    theme_classic() +
    facet_grid(. ~ method + cell_type, drop = FALSE) +
    theme(axis.text.y=element_text(size=14,color="black"), 
          axis.text.x = element_text(size = 18, colour = "black"), 
          legend.title = element_text(size = 18, colour = "black"), 
          axis.title = element_text(size = 18, colour = "black"),
          legend.text = element_text(size = 18, colour = "black"),
          axis.title.y = element_blank())

combined_figure.padj
ggsave("combined.Reactome.padj.nOverlap.pdf", combined_figure.padj, width = 18, height = 7, useDingbats = FALSE, path = "./Paper_figures/")
ggsave("combined.Reactome.padj.nOverlap.png", combined_figure.padj, width = 18, height = 7, device = "png", path = "./Paper_figures/")

