#$ -cwd
#$ -o ../../processed-data/03_SPEAQeasy_Y2/01-make_manifest.log
#$ -e ../../processed-data/03_SPEAQeasy_Y2/01-make_manifest.log
#$ -N make_manifest

module load conda_R/4.1.x
Rscript 01-make_manifest.R
