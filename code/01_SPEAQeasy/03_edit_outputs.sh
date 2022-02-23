#$ -cwd
#$ -o ../../processed-data/01_SPEAQeasy/03_edit_outputs.log
#$ -e ../../processed-data/01_SPEAQeasy/03_edit_outputs.log
#$ -l mf=20G,h_vmem=20G
#$ -N edit_outputs_03

module load conda_R/4.1.x
Rscript 03_edit_outputs.R
