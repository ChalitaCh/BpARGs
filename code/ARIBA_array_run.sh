#!/bin/bash

#PBS -l walltime=05:00:00
#PBS -l select=1:ncpus=2:mem=4gb
#PBS -J 1-3341


#Start of the script

cd $PBS_O_WORKDIR
source ~/anaconda3/bin/activate

conda activate BP

#Get the array index id

N=$(sed -n "${PBS_ARRAY_INDEX}p" samplenames_3341.txt)

#Select the output directory name and ARIBA database name
OUT="OUTPUT_RES"
DATA="ref_res"

#Running ariba

echo "Running ariba"

ariba run --threads 2 ARIBA/$DATA ../../projects/chalitaproject/live/RAWREADS/${N}_1.fastq.gz ../../projects/chalitaproject/live/RAWREADS/${N}_2.fastq.gz $HOME/$OUT/$N

echo "Finish running, exiting with no problem?"

#End of the script
