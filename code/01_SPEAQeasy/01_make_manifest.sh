#$ -cwd
#$ -o ../../processed-data/01_SPEAQeasy/01_make_manifest.log
#$ -e ../../processed-data/01_SPEAQeasy/01_make_manifest.log
#$ -N make_manifest_01

module load conda_R/4.1.x
Rscript 01_make_manifest.R
