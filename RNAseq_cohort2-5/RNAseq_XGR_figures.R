#Aim to draw graph for XGR output

library(tidyverse)
library(ggrepel)

#Define colours
Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000") 
#CD4 009E73 Okabe_Ito[3]
#CD8 0072B2 Okabe_Ito[5]
#CD14 D55E00 Okabe_Ito[6]

#Read in files
RNA.CD4.Reactome <- read.table("./Reactome/RNA.CD4.sig.Reactome.txt",
                               header = TRUE, sep = "\t")
RNA.CD8.Reactome <- read.table("./Reactome/RNA.CD8.sig.Reactome.txt",
                               header = TRUE, sep = "\t")
RNA.CD14.Reactome <- read.table("./Reactome/RNA.CD14.sig.Reactome.txt",
                               header = TRUE, sep = "\t")

RNA.CD4.KEGG <- read.table("./KEGG/RNA.CD4.sig.KEGG.txt",
                               header = TRUE, sep = "\t")
RNA.CD8.KEGG <- read.table("./KEGG/RNA.CD8.sig.KEGG.txt",
                               header = TRUE, sep = "\t")
RNA.CD14.KEGG <- read.table("./KEGG/RNA.CD14.sig.KEGG.txt",
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

RNA.CD4.Reactome.select <- make_df1(RNA.CD4.Reactome, "CD4")
RNA.CD8.Reactome.select <- make_df1(RNA.CD8.Reactome, "CD8")
RNA.CD14.Reactome.select <- make_df1(RNA.CD14.Reactome, "CD14")
RNA.CD4.KEGG.select <- make_df1(RNA.CD4.KEGG, "CD4")
RNA.CD8.KEGG.select <- make_df1(RNA.CD8.KEGG, "CD8")
RNA.CD14.KEGG.select <- make_df1(RNA.CD14.KEGG, "CD14")

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

RNA.Reactome <- join_fn(RNA.CD4.Reactome.select, RNA.CD8.Reactome.select, RNA.CD14.Reactome.select)
RNA.KEGG <- join_fn(RNA.CD4.KEGG.select, RNA.CD8.KEGG.select, RNA.CD14.KEGG.select)



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

RNA.CD4.Reactome.select2 <- make_df2(RNA.CD4.Reactome, "CD4")
RNA.CD8.Reactome.select2 <- make_df2(RNA.CD8.Reactome, "CD8")
RNA.CD14.Reactome.select2 <- make_df2(RNA.CD14.Reactome, "CD14")
RNA.CD4.KEGG.select2 <- make_df2(RNA.CD4.KEGG, "CD4")
RNA.CD8.KEGG.select2 <- make_df2(RNA.CD8.KEGG, "CD8")
RNA.CD14.KEGG.select2 <- make_df2(RNA.CD14.KEGG, "CD14")

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

RNA.Reactome2 <- join_factor_fn(RNA.CD4.Reactome.select2, RNA.CD8.Reactome.select2, 
                       RNA.CD14.Reactome.select2, RNA.Reactome)
RNA.KEGG2 <- join_factor_fn(RNA.CD4.KEGG.select2, RNA.CD8.KEGG.select2, 
                                RNA.CD14.KEGG.select2, RNA.KEGG)

#Draw plots of name vs log2fc with size as -log10padj
plot_fn <- function (df, scale){
    plot <- ggplot(df, aes(x = name, y = log2fc, size = minuslog10padj))+
        geom_point(aes(colour = cell_type)) +
        scale_size(range = c(1, 12))+
        scale_y_continuous(name = "log2(FC)", limits = scale) +
        scale_colour_manual(values = c(CD4 = Okabe_Ito[3], CD8 = Okabe_Ito[5], CD14 = Okabe_Ito[6])) +
        coord_flip() +
        theme_classic() + 
        facet_wrap(~cell_type) + 
        theme(axis.text.y=element_text(size=24,color="black"), 
              axis.text.x = element_text(size = 24, colour = "black"), 
              legend.title = element_text(size = 24, colour = "black"), 
              axis.title = element_text(size = 24, colour = "black"),
              legend.text = element_text(size = 24, colour = "black"),
              axis.title.y = element_blank(),
              strip.text = element_blank())
        return (plot)
}

Reactome.plot.padj <- plot_fn(RNA.Reactome2, scale = c(0, 4))
KEGG.plot.padj <- plot_fn(RNA.KEGG2, scale = c(0, 5))
Reactome.plot.padj
KEGG.plot.padj

ggsave("RNA.Reactome.padj.pdf", Reactome.plot.padj, width = 20, height = 12, useDingbats = FALSE, path = "./figures/")
ggsave("RNA.Reactome.padj.png", Reactome.plot.padj, width = 20, height = 12, device = "png", path = "./figures/")
ggsave("RNA.KEGG.padj.pdf", KEGG.plot.padj, width = 20, height = 12, useDingbats = FALSE, path = "./figures/")
ggsave("RNA.KEGG.padj.png", KEGG.plot.padj, width = 20, height = 12, device = "png", path = "./figures/")


#Draw the plot with y axis log2f and size -nOverlap
Reactome.plot.nOverlap <- ggplot(RNA.Reactome2, aes(x = name, y = log2fc, size = nOverlap)) + 
    geom_point(aes(colour = cell_type)) +
    scale_y_continuous(name = "log2(FC)", limits = c(0, 3.5)) +
    scale_colour_manual(values = c(CD4 = Okabe_Ito[3], CD8 = Okabe_Ito[5], CD14 = Okabe_Ito[6])) +
    coord_flip() +
    theme_classic() + 
    facet_wrap(~cell_type) + 
    theme(axis.text.y=element_text(size=24,color="black"), 
          axis.text.x = element_text(size = 24, colour = "black"), 
          legend.title = element_text(size = 24, colour = "black"), 
          axis.title = element_text(size = 24, colour = "black"),
          legend.text = element_text(size = 24, colour = "black"),
          axis.title.y = element_blank(),
          strip.text = element_blank())

Reactome.plot.nOverlap

ggsave("RNA.Reactome.nOverlap.pdf", Reactome.plot.nOverlap, width = 12, height = 14, useDingbats = FALSE, path = "./figures/")
ggsave("RNA.Reactome.nOverlap.png", Reactome.plot.nOverlap, width = 12, height = 14, device = "png", path = "./figures/")

#Draw a plot of FDR vs z-score with size by number of genes


Reactome.plot.zscore <- ggplot(RNA.Reactome2, aes(x = zscore, y = minuslog10padj, size = nOverlap)) + 
    geom_point(aes(colour = cell_type)) +
    scale_colour_manual(values = c(CD4 = Okabe_Ito[3], CD8 = Okabe_Ito[5], CD14 = Okabe_Ito[6])) +
    scale_y_continuous(name = "-log10(padj)", limits = c(0, 10)) +
    scale_x_continuous(name = "Z-score", limits = c(0, 10)) +
    theme_classic() +
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

ggsave("RNA.Reactome.zscore.pdf", Reactome.plot.zscore, width = 12, height = 10, useDingbats = FALSE, path = "./figures/")
ggsave("RNA.Reactome.zscore.png", Reactome.plot.zscore, width = 12, height = 10, device = "png", path = "./figures/")

#Draw a plot of FDR vs z-score with size by number of genes for KEGG


KEGG.plot.zscore <- ggplot(RNA.KEGG2, aes(x = zscore, y = minuslog10padj, size = nOverlap)) + 
    geom_point(aes(colour = cell_type)) +
    scale_colour_manual(values = c(CD4 = Okabe_Ito[3], CD8 = Okabe_Ito[5], CD14 = Okabe_Ito[6])) +
    scale_y_continuous(name = "-log10(padj)", limits = c(2,8)) +
    scale_x_continuous(name = "Z-score", limits = c(0, 12)) +
    theme_classic() +
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
KEGG.plot.zscore

ggsave("RNA.KEGG.zscore.pdf", KEGG.plot.zscore, width = 12, height = 10, useDingbats = FALSE, path = "./figures/")
ggsave("RNA.KEGG.zscore.png", KEGG.plot.zscore, width = 12, height = 10, device = "png", path = "./figures/")
