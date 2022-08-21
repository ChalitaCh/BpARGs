#!/bin/bash

#PBS -l walltime=05:00:00
#PBS -l select=1:ncpus=1:mem=4gb
#PBS -J 1-3341

#Start of the script

cd $PBS_O_WORKDIR

#Get the array index id

N=$(sed -n "${PBS_ARRAY_INDEX}p" samplenames_3341.txt)

#Starting of the script

#Copy the outputs from ARIBA to another directory
cp OUTPUT_PUT/${N}/assembled_genes.fa.gz GENE_ASSEMBLED_PUT/${N}_genes.fa.gz

#Remove the unnecessary part of the sequence name (only the gene name left)
zcat GENE_ASSEMBLED_PUT/${N}_genes.fa.gz | sed 's/[.].*$//' > GENE_ASSEMBLED_PUT_EDITED/${N}_genes_edited.fa

#gzip the file to safe the storage space 
gzip GENE_ASSEMBLED_PUT_EDITED/${N}_genes_edited.fa

#make directory for each sample

mkdir -p GENE_ASSEMBLED_EDITED/${N}

cd GENE_ASSEMBLED_EDITED/${N}

#splite the fasta files into the individual genes
awk '/^>/ {out = substr($1, 2) ".fasta"; print > out} !/^>/ {print >> out}' <(gzip -dc ../${N}_genes_edited.fa.gz)

#End of the script