#Characterised genes SNPs exploration and visualisation

#remove previous session variables

rm(list = ls())

#Load require packages

require(ggplot2)
require(tidyverse)
require(stringr)


#Load datasets and metadata

#Metadata
metadata <- read.csv("../data/BP.3341.metadata.source.edited2.csv", header = TRUE)

#Gene information
gene_info <- read.csv("../data/BP_ARGs_info.csv", header = TRUE)

#Ariba cluster name list

ariba_cluster <- read.delim("../data/Ariba_ref_clusters.tsv", header = FALSE)

ariba_cluster$V2 <- c("bpeB", "bpeF","bpeA","bpeE","amrB", "bpeS", "amrA", "omp38",
                  "bpeT", "BPSL3085", "folP", "penA", "pbp3_1", "blaTEM",
                  "blaOXA", "oxa_57","pbp1_2", "ptr1", "lrlR", "amrR",
                  "bpeR","BCAS0292", "folA","pbp3_2", "BPSS0262",
                  "II2144", "II2141", "gyrA", "pbp1_5", "pbp1_4",
                  "pbp2", "pbp1_3", "pbp1_1", "pbp3_3", "oprA", "oprB",
                  "oprC", "norM", "pbp6")

#all snps
all_snps <- read.delim("../results/Ariba_all_results.txt", header = TRUE)

#Add the associated metadata

all_snps$Country <- metadata$Country[match(all_snps$sample, metadata$name)]
all_snps$Continent <- metadata$Continent[match(all_snps$sample, metadata$name)]
all_snps$gene_name <- ariba_cluster$V2[match(all_snps$cluster, ariba_cluster$V1)]

all_snps$ref_ctg_effect <- if_else(all_snps$ref_ctg_change == "." & all_snps$ref_ctg_effect == ".", "NO_CHANGE", all_snps$ref_ctg_effect)
all_snps$ref_ctg_effect <- if_else(all_snps$ref_ctg_effect == ".", "SYN", all_snps$ref_ctg_effect)

#Calculate the unique number of variation types observed in each characterised gene

var_freq <- all_snps %>%
  group_by(gene_name, ref_ctg_effect) %>%
  summarise(n = length(unique(ref_ctg_change)))

#make it into a wide format for visualisation in the report
var_freq_wide <- spread(var_freq, ref_ctg_effect, n)
#change NA values to 0
var_freq_wide[is.na(var_freq_wide)] <- 0

write.csv("../results/BP_characterised_variants_type.csv", quote = FALSE, row.names = FALSE)


#Plotting the gene coordiantion and SNPs distribution on the gene sequence
gene_list <- unique(all_snps$cluster)

prob_idx <- c()

for (i in 1:length(gene_list)) {
  
  tryCatch(
    {
      
      graphics.off()
      cluster_index <- which(ariba_cluster$V1 == gene_list[i])
      gene_index <- which(gene_info$Other.name == ariba_cluster$V2[cluster_index])
      gene_start <- gene_info$Start[gene_index]
      gene_end <- gene_info$End[gene_index]
      gene_name <- gene_info$Other.name[gene_index]
      
      data_syn <- all_snps %>%
        filter(cluster == gene_list[i], ref_ctg_effect == "SYN" ) %>%
        group_by(ref_ctg_change, ref_start) %>%
        summarise(Percent = (n()/3341)*100)
      
      data_syn$start_nt <- as.numeric(data_syn$ref_start) + gene_start
      
      data_nonsyn <- all_snps %>%
        filter(cluster == gene_list[i], ref_ctg_effect == "NONSYN") %>%
        group_by(ref_ctg_change,ref_start) %>%
        summarise(Percent = (n()/3341)*100)
      
      data_nonsyn$start_nt <- as.numeric(data_nonsyn$ref_start) + gene_start

      data_to_label <- data_nonsyn %>%
        filter(Percent > 2)
      
      pdf(file = paste("../results/SNPS/", gene_name, "_snps.pdf"),
          width = 8.27,
          height = 5.83)
      
      plot(c(gene_start - 1000, gene_end + 1000),
           c(0,110), 
           type = "n",
           xlab = gene_name , 
           ylab = "",
           bty='n', yaxt="n")
      
      rect(xleft = gene_start, 
           ybottom = 0,
           xright = gene_end, 
           ytop = 10, 
           border = "black")
      
      segments(x0 = data_syn$start_nt,
               y0 = 0, 
               x1 = data_syn$start_nt, 
               y1= data_syn$Percent,
               col = "blue")
      segments(x0 = data_nonsyn$start_nt, 
               y0 = 0, 
               x1 = data_nonsyn$start_nt, 
               y1= data_nonsyn$Percent,
               col = "red")
      text(x = data_to_label$start_nt,
           y = data_to_label$Percent + 5,
           label = data_to_label$ref_ctg_change)
      
      dev.off()
    }, error = function(e){
      prob_idx <- c(prob_idx, i)
    }
  )
  
}
