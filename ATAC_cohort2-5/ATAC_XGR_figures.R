#Aim to draw graph for XGR output
setwd("~/R/XGR/ATAC_cohort2-5")
library(tidyverse)
library(ggrepel)

#Define colours
Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000") 
#CD4 009E73 Okabe_Ito[3] blueish green
#CD8 0072B2 Okabe_Ito[5] blue
#CD14 D55E00 Okabe_Ito[6] vermillion

#Read in files
ATAC.CD4.Reactome <- read.table("./Reactome/ATAC.CD4.Reactome_2.txt",
                               header = TRUE, sep = "\t")
ATAC.CD8.Reactome <- read.table("./Reactome/ATAC.CD8.Reactome_2.txt",
                               header = TRUE, sep = "\t")
ATAC.CD14.Reactome <- read.table("./Reactome/ATAC.CD14.Reactome_2.txt",
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

#Draw plots of name vs log2fc with size as -log10padj
plot_fn <- function (df, scale){
    plot <- ggplot(df, aes(x = name, y = log2fc, size = minuslog10padj))+
        geom_point(aes(colour = cell_type)) +
        scale_y_continuous(name = "log2(FC)", limits = scale) +
        scale_colour_manual(values = c(CD4 = Okabe_Ito[3], CD8 = Okabe_Ito[5], CD14 = Okabe_Ito[6])) +
        coord_flip() +
        theme_classic() + 
        facet_wrap(~cell_type, drop = FALSE) + 
        theme(axis.text.y=element_text(size=24,color="black"), 
              axis.text.x = element_text(size = 24, colour = "black"), 
              legend.title = element_text(size = 24, colour = "black"), 
              axis.title = element_text(size = 24, colour = "black"),
              legend.text = element_text(size = 24, colour = "black"),
              axis.title.y = element_blank(),
              strip.text = element_blank())
    return (plot)
}

Reactome.plot.padj <- plot_fn(ATAC.Reactome2, scale = c(0, 6))
Reactome.plot.padj

ggsave("ATAC.Reactome.padj.pdf", Reactome.plot.padj, width = 15, height = 6, useDingbats = FALSE, path = "./figures/")
ggsave("ATAC.Reactome.padj.png", Reactome.plot.padj, width = 15, height = 6, device = "png", path = "./figures/")
