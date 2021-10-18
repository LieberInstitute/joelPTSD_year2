#  As of 10/18/21, identify which original FASTQ files look suspicious, and
#  check if the "recopy" provided "correct" versions of the corrupted files.
#  This script actually finds that the recopy was only partially complete and
#  does not cover all corrupted original files.

library('jaffelab')
library('here')

fastq_src_orig = '/dcl01/ajaffe/data/Nina/Joel_R01/fastq/Year2'
fastq_src_new = '/dcl01/ajaffe/data/Nina/Joel_R01/fastq/Year2/HCHY3BBXY_recopy'

#  The files that look suspicious so far are those that are very small in size
suspicious_files = basename(
    system(
        paste0('ls -s ', fastq_src_orig, "/*.fastq.gz | awk '$1 < 300000 {print $2}'"),
        intern = TRUE
    )
)

#  Compare md5 sums for suspicious vs. recopied files
md5_list = list()
for (f in suspicious_files) {
    stopifnot(file.exists(file.path(fastq_src_new, f)))
    
    md5_list[[f]] = system(paste('md5sum', file.path(fastq_src_orig, f), file.path(fastq_src_new, f)), intern = TRUE)
}
