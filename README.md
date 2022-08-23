#  Antibiotic resistance genes analyses in *Burkholderia pseudomallei*

This repository is for the analyses of Antibiotic Resistance Genes (ARGs) evolutionary dynamics in *B. pseudomallei* using a global collection of the bacteria 3,341 whole-genome sequences from the public databases. This project is part of the MSc Computational methods in ecology and evolution at Imperial College London.

## Description

This project aims to identify both novel and characterised ARGs, and the prevalence of resistant variants in a global population of *B. pseudomallei*. Furthermore, I hope to disentangle how the resistance genes and variants bring about, and their association to human activities such as antibiotic treatment in clinical settings. This project utilised the combination of R and shell script in the analyses. The directories are listed as; data, code, results, and sandbox. 

## Getting Started

### Dependencies
All the softwares and packages used in this project are listed. The dependencies for installation are according to the recommended guideline of each software.
* Softwares:
    * [Unicycler](https://github.com/rrwick/Unicycler)
    * [Quast](https://github.com/ablab/quast)
    * [Prokka](https://github.com/tseemann/prokka)
    * [PopPUNK](https://poppunk.readthedocs.io/en/latest/)
    * [Mafft](https://mafft.cbrc.jp/alignment/software/)
    * [SNP-sites](https://github.com/sanger-pathogens/snp-sites)
    * [IQ-TREE](http://www.iqtree.org)
    * [ARIBA](https://github.com/sanger-pathogens/ariba)
    * [Panaroo] (https://gtonkinhill.github.io/panaroo/#/)
    * [Gubbins](https://github.com/nickjcroucher/gubbins)
    * [PAML](https://github.com/abacus-gene/paml)
    * [BactDating](https://github.com/xavierdidelot/BactDating)
* R version 4.2.0
* R packages: ggplot2, tidyverse, ape, phytools, reshape2, R.utils, ggpubr, treeio, ggtree, chopper (fas2phy), doParallel

### Code

Please ensure that the path directory is in the correct one. Some scripts might require the results from previous script before.

* ARIBA_array_run.sh : to identify ARGs in the whole-genome sequencing using ARIBA and customised database
* ARIBA_pull_genes.sh : to pull the ARGs identified from ARIBA result for further analyses, creating a single fasta file of the ARGs.
* ARIBA_summary.sh : to contanate all the genetic variations results from ARIBA together for further analyses
* codeml.ctl : an example of control file used to test for selection that acting on the gene using PAML. To run this analysis, please ensure that the sequence alignment, gene phylogeny and the control file are all in the same directory before execute the following command;
```
../paml4.8/bin/codeml
```
* genes_process.sh : to create a gene sequence alignment and phylogeny generated from `ARIBA_pull_genes.sh` and later used in conjunction with `codeml.ctl`
* Gubbins_process.sh : to generate the frequency table of recombination events happended in each genome coordinate
* Gubbins_visualise.R : to identify the recombination hotspots and genes present in them.
* Gubbins_enrichment.R : to test for significant likelihood of the ARGs and mobile genetic elements present in the recombination hotspots compared to non-hotspots in each lineage.
* mole_dating_function.R : a divergent-time estimation function, developing from the `bactdate` function in BactDating package to run using multicore cpu.
* mole_dating_parallel.R : to run the molecular dating using multicore cpu.
* Overall_dataset_F1B.R : R script to plot the figure 1B in the report, representing the data distribution across year and sources.
* PAML_LK : to calculate the log-likelihood ratio statistic test and AIC score of the models'log-likelihood.
* SNPs_visualisation_FS1_TS3.R : R script to visualise the dirtribution of SNPs variantions seen in the dataset, both along the nucleotide sequence and the unique number of variants found (Figure S1 and Table S3). To run this script, please identify the variants of each isolate using `Ariba_array_run.sh`, then merge the results together using `Ariba_summary.sh`. The file name and absolute path may need to be changed accordingly.
* Time_trees_F3.R : Time-calibrated phylogeny's visualisation with their associated data/variants (Figure 3).


## Author

Chalita Chomkatekaew
[@CChomkatekaew](https://twitter.com/CChomkatekaew)
