library("jaffelab")
library("here")
library("sessioninfo")

#  Paths to input and output FASTQs, and output manifest
fastq_src <- "/dcs04/lieber/data/transfers/Joel_R01_Y2_fourthtry/FASTQ"
fastq_dest <- here("raw-data", "FASTQ_Y2")
man_path <- here("processed-data", "03_SPEAQeasy_Y2", "samples.manifest")

#  Create directories that don't exist
dir.create(fastq_dest, showWarnings = FALSE)

#  The original paths to the first and second reads
r1_src <- list.files(
    fastq_src, ".*_L00._R1_001\\.fastq\\.gz", full.names = TRUE
)
r2_src <- list.files(
    fastq_src, ".*_L00._R2_001\\.fastq\\.gz", full.names = TRUE
)
stopifnot(length(r1_src) == length(r2_src))

ids <- ss(basename(r1_src), "_")

#  The destination paths (where we will place symbolic links)
r1_dest <- file.path(fastq_dest, basename(r1_src))
r2_dest <- file.path(fastq_dest, basename(r2_src))

#  Form the symbolic links
all(file.symlink(r1_src, r1_dest))
all(file.symlink(r2_src, r2_dest))

#  Create and write the manifest
man <- paste(r1_dest, 0, r2_dest, 0, ids, sep = "\t")
writeLines(man, con = man_path)

session_info()
