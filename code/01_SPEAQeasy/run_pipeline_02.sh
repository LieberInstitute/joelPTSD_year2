#!/bin/bash
#$ -l bluejay,mem_free=40G,h_vmem=40G,h_fsize=800G
#$ -o ../../processed-data/01_SPEAQeasy/SPEAQeasy_output.log
#$ -e ../../processed-data/01_SPEAQeasy/SPEAQeasy_output.log
#$ -cwd

ORIG_DIR=/dcl01/lieber/ajaffe/Nick/SPEAQeasy

#  Get absolute path to 'joelPTSD_year2' repo
base_dir=$(git rev-parse --show-toplevel)

module load nextflow
export _JAVA_OPTIONS="-Xms8g -Xmx10g"

nextflow $ORIG_DIR/main.nf \
    --sample "paired" \
    --reference "hg38" \
    --strand "reverse" \
    --annotation "/dcl01/lieber/ajaffe/Nick/SPEAQeasy/Annotation" \
    -with-report "${base_dir}/processed-data/01_SPEAQeasy/execution_reports/JHPCE_run.html" \
    -w "${base_dir}/processed-data/01_SPEAQeasy/work" \
    --input "${base_dir}/processed-data/01_SPEAQeasy" \
    --output "${base_dir}/processed-data/01_SPEAQeasy/pipeline_output" \
    -profile ptsd_y2

#  Produces a report for each sample tracing the pipeline steps
#  performed (can be helpful for debugging).
#
#  Note that the reports are generated from the output log produced in the above
#  section, and so if you rename the log, you must also pass replace the filename
#  in the bash call below.
echo "Generating per-sample logs for debugging..."
bash $ORIG_DIR/scripts/generate_logs.sh ../../processed-data/01_SPEAQeasy/SPEAQeasy_output.log
