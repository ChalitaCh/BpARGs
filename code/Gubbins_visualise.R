#Recombination analyses
#This script is for the visualisation of the recombination frequency across the chromosomes of B. pseudomallei
#The plot is then annotated with the known and characteriesed ARGs and MGEs (known location on the K96243 reference genome)

#remove any object in R before starting

rm(list = ls())
graphics.off()

#packages required for analyses

require(ggplot2)
require(tidyverse)
require(R.utils)
require(ggpubr)

#load the dataset(s) required

    #All the genes annotation in Bp K96243 reference genome
    BP_genes <- read.csv("../data/ref_BP_genes_annotation_clean.csv", header = TRUE)

    #Finding all the output files from the Gubbins_process.sh

    main_dir <- "../results/GUBBINS/CHROM1/"

    files_list <- list.files(path = main_dir, pattern = "\\.recom.txt$", recursive = TRUE )

    #Specificied the lineages for plotting the figure
    #files_list <- files_list[c(14,18,22,2,6,9,15)]

    #Select genes for specific chromosome, please change this accordingly to what chromosome you are working with
        #Chr1 : NC_006350.1
        #     : (0,4074542,200000)
        #Chr2 : NC_006351.1
        #     : (0,3100174,200000)

    all_chr <- BP_genes %>%
        filter(seq_id == "NC_006350.1"  ) 
  
    #Create empty dataframes to store the results

    all_genes_found_chr1 <- data.frame(GeneID.K96243 = as.character(),
                              Other.name = as.character(),
                              Strand = as.character(),
                              Start = as.numeric(),
                              End = as.numeric(),
                              Annotate = as.character(),
                              Function = as.character(),
                              Cluster = as.character(),
                              No_regions = as.numeric()
    )

    all_genes_found_chr2 <- data.frame(GeneID.K96243 = as.character(),
                              Other.name = as.character(),
                              Strand = as.character(),
                              Start = as.numeric(),
                              End = as.numeric(),
                              Annotate = as.character(),
                              Function = as.character(),
                              Cluster = as.character(),
                              No_regions = as.numeric()
    )
  
    plot_chr1 <- list()
    plot_chr2 <- list()
  
    cluster_name <- c()
  
#Find the genome coordination with recombination frequency higher than the 4th quantile
   
for (i in 1:length(files_list)) {

  #get the lineage name from the file name
  name <- unlist(strsplit(files_list[i], split = "[.]+"))
  
  #save the lineage name together
  cluster_name <- c(cluster_name, name[1])
  
  #read the recombination frequency data
  clust <- read.delim(paste(main_dir,files_list[i], sep = ""), header = FALSE)
  
  #identify the 95th quantile of recombination events in that lineage
  q95 <- quantile(clust$V2, c(0.95))

  #identify the location of which it has higher recombination frequency to the 95th percentile
  seq_pos <- clust$V1[which(clust$V2 >= q95)]

  #Change the location vectors into the interval dataframe
  recombine_regions <- as.data.frame((seqToIntervals(sort(seq_pos))))
  
  #Expand the recombination regions to cover 10kb neighboring regions of the hotspots
  recombine_regions_wide <- recombine_regions
  recombine_regions_wide$from <- recombine_regions_wide$from - 10000
  recombine_regions_wide$to <- recombine_regions_wide$to + 10000
  
  #All genes in BP ref annotation that located with in the recombination hotspots
  
  all_pos <- intersect(all_chr$start, unlist(apply(recombine_regions_wide, 1, function(a) a[1]:a[2])))
  
  all_name <- subset(all_chr, start %in% all_pos) 
  
  all_name_new <- all_name %>%
    select(old_locus_tag,locus_tag,strand,start,end,product, Function) %>%
    rename(GeneID.K96243 = old_locus_tag,
           Other.name = locus_tag,
           Strand = strand,
           Start = start,
           End = end,
           Annotate = product,
           Function = Function) %>%
    mutate(Cluster = name[1],
           No_regions = no_regions)
  
  #combine all the genes found in the recombination together into a big dataframe
  all_genes_found_chr1 <- rbind(all_genes_found_chr1, all_name_new)
  
  #Plot the recombination frequency along with the MGEs + ARGs position on the chr1
  recombine_plot <- clust %>%
    ggplot(aes(x=V1,V2))+
    geom_line() +
    xlab(NULL) +
    ylab(NULL) +
    scale_x_continuous(breaks =seq(0,4074542,200000)) + #change the x-axis scale when change the chromosome
    geom_hline(yintercept = q95, color = "red") +
    annotate("rect", xmin = recombine_regions$from, xmax = recombine_regions$to, ymin = -Inf,
           ymax = Inf, alpha = 0.6, fill = "grey") +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          plot.background = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank())
  
  #annotate the ARGs and MGEs on the plot
  final_plot <- recombine_plot +
      geom_vline(aes(xintercept = Start),data = all_name_new[all_name_new$Function == "MGE",], color = "green") + #for the MGEs
      geom_vline(aes(xintercept = Start), data = all_name_new[all_name_new$Function == "ARG",], color = "blue")  #for the ARGs


  plot_chr1[[i]] <- final_plot
    
}

#save the genes found in the recombination hotspots as a csv file for further analayses

write.csv(all_genes_found_chr1, file = "../results/genes.recomb.chr1.csv", quote = FALSE, row.names = FALSE)
  
#Plotting the figures for report

#Combinding the lineages of interest into the same plot  
  plot_chr1_all <- ggarrange(plot_chr1[[2]],
                        plot_chr1[[4]],
                        plot_chr1[[5]],
                        plot_chr1[[6]],
                        plot_chr1[[7]],
                        nrow = 5, ncol = 1,
                        labels = c("Bp7","Bp12",
                                   "Bp16","Bp19","Bp21"))
  
  plot_chr2_all <- ggarrange(plot_chr2[[2]],
                        plot_chr2[[4]],
                        plot_chr2[[5]],
                        plot_chr2[[6]],
                        plot_chr2[[7]],
                        nrow = 5, ncol = 1)

#Merge both chromosome plots together before saving
 
  plot_merge <- ggarrange(plot_chr1_all,
                        plot_chr2_all,
                        nrow = 1, ncol = 2,
                        align = "h")

#Save the plots

ggsave(plot_merge, file = "../results/all_recomb_plot.pdf" ,
       width = 16.5, height = 11.75, units = "in" )

