#$ -cwd
#$ -o ../../processed-data/01_SPEAQeasy/04_combine_years.log
#$ -e ../../processed-data/01_SPEAQeasy/04_combine_years.log
#$ -l bluejay,mf=10G,h_vmem=10G,h_fsize=50G
#$ -N combine_years_04

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load conda_R/4.1.x
Rscript 04_combine_years.R

echo "**** Job ends ****"
date
