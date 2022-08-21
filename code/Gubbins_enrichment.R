#Recombination enrichment analyses

#This script is for the identification of which lineages are significantly enriched in the ARGs and MGEs at the recombination hotspots.


#remove any object in R before starting

rm(list = ls())
graphics.off()

#packages required for analyses

require(ggplot2)
require(tidyverse)

#load the dataset(s) required

BP_genes <- read.csv("../data/ref_BP_genes_annotation_clean.csv", header = TRUE)
all_recomb_chr1 <- read.csv("../results/genes.recomb.chr1.csv", header = TRUE)
all_recomb_chr2 <- read.csv("../results/genes.recomb.chr2.csv", header = TRUE)

#Combine the genes identified in the two chromosomes together

all_recomb <- rbind(all_recomb_chr1, all_recomb_chr2)

#Calculate the probability of success in the single trial i.e. probability to find the ARGs/MGEs in the whole genome

    #No of genes in the whole genomes
    all_genes_no <- length(unique(BP_genes$locus_tag))

    #No of ARGs and MGEs in the whole genomes
    ARG_no <- nrow(BP_genes[BP_genes$Function %in% c("MGE","ARG"),])

    #the chance of getting the ARG + MGE
    p <- ARG_no/all_genes_no

#Find the likelihood of finding the ARGs and MGEs across the genome

    binom_genome <- dbinom(ARG_no, size = all_genes_no, prob = p, log = TRUE)

    conf_genome <- prop.test(x = ARG_no, n = all_genes_no, conf.level = 0.95, correct = FALSE) 

#Find the likelihood of finding the ARGs and MGEs in the recombinaiton hotspots in each lineage

cluster <- unique(all_recomb$Cluster)

#Initialise the data frame to store the results
results <- data.frame(cluster = as.character(),
                      no_ARG_hotspot = as.numeric(),
                      no_genes_hotspot = as.numeric(),
                      prob = as.numeric(),
                      conf_level_low = as.numeric(),
                      conf_level_high = as.numeric(),
                      p_value = as.numeric()
                      )

for (i in 1:length(cluster)) {
  
  #Select only the genes found in the hotspot of one lineage
  data_subset <- subset(all_recomb, Cluster == cluster[i])
  
  no_ARG_hotspot <- nrow(data_subset[data_subset$Function %in% c("MGE", "ARG"),])
  
  no_genes_hotspot <- nrow(data_subset)

  #Find the likelihood of finding ARGs and MGEs in the hostspot.
  binom_test <- dbinom(no_ARG_hotspot, size = no_genes_hotspot, prob = p, log = TRUE)
  
  conf_test <- prop.test(x = no_ARG_hotspot, n = no_genes_hotspot, conf.level = 0.95, correct = FALSE)
  
  #save all the results in the temp data frame
  df_temp <- data.frame(cluster = cluster[i],
                        no_ARG_hotspot = no_ARG_hotspot,
                        no_genes_hotspot = no_genes_hotspot,
                        prob = binom_test,
                        conf_level_low = conf_test$conf.int[1],
                        conf_level_high = conf_test$conf.int[2],
                        p_value = conf_test$p.value
                        )
  #save the results together in the final database
  results <- rbind(results, df_temp)
}

#Add the results from the whole genome analysed earlier
results[nrow(results) + 1,] <- c("Genome",
                                 ARG_no, 
                                 all_genes_no,
                                 binom_genome,
                                 conf_genome$conf.int[1],
                                 conf_genome$conf.int[2],
                                 conf_genome$p.value)

#Make sure that all are numeric before the comparison
results$prob <- as.numeric(results$prob)
results$no_ARG_hotspot <- as.numeric(results$no_ARG_hotspot)
results$no_genes_hotspot <- as.numeric(results$no_genes_hotspot)
results$conf_level_low <- as.numeric(results$conf_level_low)
results$conf_level_high <- as.numeric(results$conf_level_high)

#Compare the likelihood of finding the ARGs and MGEs between the non-hotspots and hotspots
results$outcome <- ifelse(results$prob > binom_genome, "enriched", "not")

#Plotting the results
ggplot(results, aes( x= cluster, y =prob, color = outcome)) +
  geom_point() +
  geom_errorbar(aes( ymin = prob - conf_level_low, ymax = prob + conf_level_high)) +
  theme_bw()
  theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
  labs(x = "lineage",
       y = "log-likelihood") 

