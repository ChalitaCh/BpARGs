#R script to visualise the time-calibrated phylogeny with associated data

#clean the working environment
rm(list = ls())
graphics.off()

#load required package(s) 

require(tidyverse)
require(ggplot2)
require(treeio)
require(ggtree)
require(cowplot)

#load the dataset(s)

metadata <- read.csv("../data/BP.3341.metadata.source.edited2.csv", header = TRUE)
all_snps <- read.delim("../results/Ariba_all_results.txt", header = TRUE)
gene_info <- read.delim("../data/Ariba_ref_clusters.tsv", header = FALSE)
#data wrangling

all_snps$ref_ctg_effect <- if_else(all_snps$ref_ctg_change == "." & all_snps$ref_ctg_effect == ".", "NO_CHANGE", all_snps$ref_ctg_effect)
all_snps$ref_ctg_effect <- if_else(all_snps$ref_ctg_effect == ".", "SYN", all_snps$ref_ctg_effect)

#clean the value in the column
gene_info$V2 <- c("bpeB", "bpeF","bpeA","bpeE","amrB", "bpeS", "amrA", "omp38",
                  "bpeT", "BPSL3085", "folP", "penA", "pbp3_1", "blaTEM",
                  "blaOXA", "oxa_57","pbp1_2", "ptr1", "lrlR", "amrR",
                  "bpeR","BCAS0292", "BPSL2476","pbp3_2", "BPSS0262",
                  "II2144", "II2141", "gyrA", "pbp1_5", "pbp1_4",
                  "pbp2", "pbp1_3", "pbp1_1", "pbp3_3", "oprA", "oprB",
                  "oprC", "norM", "pbp6")

#match the ariba gene name to the commonly known gene name

all_snps$gene_name <- gene_info$V2[match(all_snps$cluster, gene_info$V1)]

#select the colour for source annotation
Source.colour <- c("human" = "#77AADD",
                   "animal" = "#EE8866",
                   "env" = "#44BB99",
                   "NA" = "#EEDD88")

#load the time-calibrated tree of cluster 2_3 or lineage 4

tree <- read.beast("../results/clust2_3.chr1.dated.tree")
load("../results/clust2_3.variants.rds")

data <- metadata %>%
  filter(cluster_poppunk_updated == "2_3")

#the followings are the steps to create a presence/absence variant matrix of each gene
data$name_Sanger <- paste0(data$name_Sanger, "_1")

# data_snps <- subset(all_snps, sample %in% data$name)
# 
# cluster_name <- sort(unique(data_snps$gene_name))
# 
# clust_arg <- list()
# 
# for (i in 1:length(cluster_name)) {
# 
#   df <- as.data.frame(with(data_snps[data_snps$gene_name == cluster_name[i],], table(sample , ref_ctg_change)) > 0L) +0L
# 
#   rownames(df) <- data$name_Sanger[match(rownames(df),data$name)]
# 
#   df[sapply(df, is.numeric)] <- lapply(df[sapply(df, is.numeric)],
#                                        as.factor)
#   df <- as.matrix(df)
# 
#   clust_arg[[i]] <- df
# 
# }
# 
# names(clust_arg) <- cluster_name
# save(clust_arg, file = "../results/clust2_3.variants.rds")

#select the variants of interest

BPSL3085_S130L <- clust_arg$BPSL3085[,4]
bpeS_V40I <- clust_arg$bpeS[,32]

to_plot <- cbind(bpeS_V40I, BPSL3085_S130L)
rownames(to_plot) <- rownames(df)

#get the most recent sampling date to plot the phylogeny
year_mrsd <- max(as.numeric(data$year), na.rm = TRUE)

#plot the tree
p <- ggtree(tree, mrsd = paste0(year_mrsd,"-01-01")) %<+% data + 
  geom_tippoint(aes(colour = Source),size=3, alpha = 0.8) +
  scale_color_manual(values = Source.colour, labels = c("Human", "Animal",
                                                        "Environment", "NA")) +
  geom_range(range='length_0.95_HPD', color='red', alpha=.6, size=2) +
  theme_tree2() +
  theme(legend.text = element_text(size = 15),
        legend.title = element_text(size = 15))

#R-square value and pvalue of temporal signal
r2 <- 0.12
pvalue <- 1e-04
Rsquare <- parse(text = paste("~R^2 ==~", r2,"~P <~", pvalue))

#plot the tree with associated data
clust2_3 <- gheatmap(p, to_plot , offset=4, width=0.8, 
         font.size=3,colnames = FALSE) +
  scale_x_ggtree(breaks = seq(1500,2018, 50)) +
  scale_fill_manual(breaks = c( "1","0"),
                    values = c( "firebrick","steelblue"),
                    name = "Presence/Absence",
                    labels = c("Presence", "Absence")) +
  theme(axis.text.x = element_text(angle = 25,
                                   vjust = 1,
                                   hjust = 1,
                                   size = 15),
        legend.position = "none") +
  geom_vline(xintercept = 2000, linetype = "dotted",lwd=1) +
  geom_vline(xintercept = 1986, linetype = "dotted",lwd=1) +
  labs(tag = Rsquare) +
  theme(plot.tag.position = c(0.1, 0.04),
        plot.tag = element_text(size = 15)) +
  #annotating the SXT and DOX timeline
  annotate("segment", x = 1975, y = 60, xend = 1999, yend = 60,
           arrow = arrow(type = "closed", length = unit(0.02, "npc"))) +
  annotate(geom="text", x=1920, y=60, label="Introduction of TMP/SXT",
           color="black", size = 6) +
  annotate("segment", x = 1930, y = 50, xend = 1985, yend = 50,
         arrow = arrow(type = "closed", length = unit(0.02, "npc"))) +
  annotate(geom="text", x=1895, y=50, label="Stop using DOX",
           color="black", size = 6) 

#save the plot
ggsave(plot = clust2_3, file = "../results/clust2_3_time_trees_heatmap.pdf",
       width = 16.5, height = 11.75, units = "in")

#plot the time-calibrated tree of cluster 17 or lineage 19
tree2 <- read.beast("../results/clust17.chr1.dated.tree")
load("../results/clust17.chr1_2.variants.rds")
data2 <- metadata %>%
  filter(cluster_poppunk_updated == "17")

data2$name_Sanger <- paste0(data2$name_Sanger, "_1")

year_mrsd2 <- max(as.numeric(data2$year), na.rm = TRUE)

p2 <- ggtree(tree2, mrsd = paste0(year_mrsd2,"-01-01")) %<+% data2 + 
  geom_tippoint(aes(colour = Source),size=3, alpha = 0.8) +
  scale_color_manual(values = Source.colour, labels = c("Human", "Animal",
                                                        "Environment", "NA")) +
  geom_range(range='length_0.95_HPD', color='red', alpha=.6, size=2) +
  theme_tree2() +
  theme(legend.text = element_text(size = 15),
        legend.title = element_text(size = 15))

#R-square and P-value
r2_2 <- 0.31
pvalue2 <- 2.91e-01
Rsquare2 <- parse(text = paste("~R^2 ==~", r2_2,"~P ==~", pvalue2))

#plot the tree
clust17 <- gheatmap(p2, var_combine[,c(2,1)] , offset=4, width=0.8, 
         font.size=3,colnames = FALSE) +
  scale_x_ggtree(breaks = seq(1980,2010, 5)) +
  scale_fill_manual(breaks = c( "1","0"),
                    values = c( "firebrick","steelblue"),
                    name = "Presence/Absence",
                    labels = c("Presence", "Absence")) +
  theme(axis.text.x = element_text(angle = 25,
                                   vjust = 1,
                                   hjust = 1,
                                   size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.position = "bottom") +
  geom_vline(xintercept = 2000, linetype = "dotted",lwd=1.5) +
  geom_vline(xintercept = 1986, linetype = "dotted",lwd=1.5) +
  labs(tag = Rsquare2) +
  theme(plot.tag.position = c(0.1, 0.06),
        plot.tag = element_text(size = 15)) +
  #annotating the SXT and DOX timeline
  annotate("segment", x = 1997.5, y = 30, xend = 2000, yend = 30,
           arrow = arrow(type = "closed", length = unit(0.02, "npc"))) +
  annotate(geom="text", x=1992, y=30, label="Introduction of TMP/SXT",
           color="black", size = 6) +
  annotate("segment", x = 1990, y = 25, xend = 1986, yend = 25,
           arrow = arrow(type = "closed", length = unit(0.02, "npc"))) +
  annotate(geom="text", x=1994, y=25, label="Stop using DOX",
           color="black", size = 6) 

#save the individual plot

ggsave(plot = clust17, file = "../results/clust17_time_trees_heatmap.pdf",
       width = 16.5, height = 11.75, units = "in")

#merge the two phylogenies together for report
plot_2 <- plot_grid( clust2_3 + theme(legend.position="none"),
                   clust17 + theme(legend.position="none"),
                   nrow = 2,
                   labels = c("(a) Lineage Bp4(chr1)", 
                   " (b) Lineage Bp19(chr1)")

)

#get a legend for one plot to make a common legend
legend_b <- get_legend(clust17 + theme(legend.position="bottom"))

#add the merged plot with the common legend.
plot_done <- plot_grid(plot_2, legend_b, ncol = 1, rel_heights = c(1, .2))

#save the plot for the report
ggsave(plot = plot_done, file = "../results/time_trees_heatmap.pdf",
       width = 11.75, height = 16.5, units = "in" )
