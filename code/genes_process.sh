#!/bin/bash

#PBS -l walltime=05:00:00
#PBS -l select=1:ncpus=4:mem=8gb
#PBS -J 5-33
#Start of the script

cd $PBS_O_WORKDIR

GENE_NAME=$(sed -n "${PBS_ARRAY_INDEX}p" gene_names.txt)

#Load the dependent environments
source ~/anaconda3/bin/activate

conda activate BP

#Start of the script

#merge all the gene fasta files of all the sample together
cat ${GENE_NAME}/*.fasta > ${GENE_NAME}.all.fasta

#Align the sequences for making a tree

mafft --auto ${GENE_NAME}.all.fasta > ${GENE_NAME}.all.aln

#echo "Successful align the genes? Manually check on seaview?"

#Create a maximum likelihood tree

iqtree -s ${GENE_NAME}.all.aln -T AUTO --threads-max 4

#End of the script