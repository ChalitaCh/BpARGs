#for running Bactdating on a single lineage using mulitcores

require(doParallel)
require(BactDating)
require(ape)
require(tidyverse)


main_dir <- "/rds/general/user/cc2320/home/BP_PROJECT/"

source(paste0(main_dir,"mole_dating_function_single.R"))

no_core <- detectCores() -1

dir_file <- "GUBBINS/CHROM1/"
cluster <- c("clust25")
data <- read.csv(paste0(main_dir,"BP.3341.metadata.source.edited2.csv"), header = TRUE)
model <- "mixedgamma"
no_iter <- 4e7
chrom <- "chr1"
output_file <- "mole_data_25_chr1_repeat5.rds"

#Running the molecular dating according to the pre-defined options above.
results <- mclapply(cluster, mole_dating, mc.cores = no_core )

#Save the results
save(results, file = paste0(main_dir,"GUBBINS/", output_file))