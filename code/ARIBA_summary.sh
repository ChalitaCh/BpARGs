#!/bin/bash

# Author: Chalita Chomkatekaew chalita.chomkatekaew20@imperial.ac.uk
# Script: Ariba_summary_all.sh
# Desc: contanate all the results together for further analyses
# Arguments: list_file.txt - name of all the samples
# Date: May 2022

LIST="samplenames_3341.txt"
RESULT_FILE_PREFIX="Ariba_all_results_putative"
OUTPUT_DIR="OUTPUT_PUT"

#For every sample; do

cat $LIST | while read name
    do

        #Fliter out only the rows with 'NONSYN' var_type tag and get the selected columns
        # 1 : ariba_ref_name
        # 7 : cluster no from Ariba output
        # 10 : %identiy between ref sequence and contigs
	# 15 : var_type
        # 19 : ref_ctg_change
        # 20 : ref_ctg_effect
        # 21 : ref_start
        # 22 : ref_end
        # 23 : ref_nt
        # 26 : ctg_nt

        awk 'NR!=1' ${OUTPUT_DIR}/${name}/report.tsv | cut -f 1,7,10,15,19,20,21,22,23,26  > sandbox/temp_output_$name.txt

        #Add the sample name in the last column of the merged file
        sed "s|$|\t$name|" "temp/temp_output_$name.txt" > sandbox/snps_final_$name.txt
    done

#remove temp files

rm temp/temp_output_*.txt

#Merge all sample files together with a new header

{ echo -e "ref_gene_name\tcluster\tpercent_id\tref_ctg_change\tref_ctg_effect\tref_start\tref_end\tref_nt\tctg_nt\tsample" ; cat sandbox/snps_final_* ; } > results/$RESULT_FILE_PREFIX.txt

#remove temp files

rm sandbox/snps_final_*.txt