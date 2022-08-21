#!/bin/bash

# Author: Chalita Chomkatekaew chalita.chomkatekaew20@imperial.ac.uk
# Script: gubbins_process.sh
# Desc: create a recombination frequency table from gubbins output
# Arguments: Gubbins output - recombination.predictions.gff file
# Date: June 2022
# The output : 1st column is genome coordinate and 2nd column is the recombination frequency

chrom="CHROM1"
cluster_name="lineage_names.txt"

cat ../results/GUBBINS/$cluster_name | while read name; do
    
    #extracting the predicted recombination regions and preparing for 'seq' command between the regions in shell script
    awk '{print "seq "$4" "$5" > seq."$4"."$5".txt"}' ../results/GUBBINS/$chrom/${name}.recombination_predictions.gff > ../results/GUBBINS/seq.${name}.sh 

    #executes the shell script
    sh ../results/GUBBINS/seq.${name}.sh

    #cat all the seq results in text files, sort them before counting the occurence/frequency
    cat seq*txt | sort | uniq -c | awk '{print $2"\t"$1}' > ../results/GUBBINS/$chrom/${name}.output.recom.txt

    rm seq*txt


done