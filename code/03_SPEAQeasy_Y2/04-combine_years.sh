#!/bin/bash
#$ -cwd
#$ -o ../../processed-data/03_SPEAQeasy_Y2/04-combine_years.log
#$ -e ../../processed-data/03_SPEAQeasy_Y2/04-combine_years.log
#$ -l mf=50G,h_vmem=50G
#$ -N combine_years

module load conda_R/4.1.x
Rscript 04-combine_years.R
