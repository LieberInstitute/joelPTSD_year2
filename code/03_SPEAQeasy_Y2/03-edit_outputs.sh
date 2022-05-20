#!/bin/bash
#$ -cwd
#$ -o ../../processed-data/03_SPEAQeasy_Y2/03-edit_outputs.log
#$ -e ../../processed-data/03_SPEAQeasy_Y2/03-edit_outputs.log
#$ -l mf=50G,h_vmem=50G
#$ -N edit_outputs

module load conda_R/4.1.x
Rscript 03-edit_outputs.R
