#   RSEs for year 1 + 2 will ultimately be merged with year 3, but this requires
#   a couple steps before the objects are compatible. Not all columns in the
#   colData are the same, year 2 uses "primary" features (a superset of the
#   "main" features used in the other years), and EntrezIDs are sometimes
#   different between objects. In all cases, we modify just the year 2 objects
#   to be compatible, then save in a new directory.

library("SummarizedExperiment")
library("here")
library("sessioninfo")

rse_dir_1and3 <- here(
    "processed-data", "01_SPEAQeasy", "updated_output", "merged"
)
rse_dir_2 <- here(
    "processed-data", "03_SPEAQeasy_Y2", "pipeline_output", "count_objects"
)

out_dir <- here("processed-data", "03_SPEAQeasy_Y2", "updated_output", "year_2")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

################################################################################
#   Functions
################################################################################

#   Convenience function for locating a particular RSE in a given directory,
#   without needing to specify the variable name of the loaded object anywhere
load_obj <- function(dir_name, var_name, verbose = TRUE) {
    rse_path <- list.files(
        dir_name,
        pattern = paste0("^", var_name, "_.*n.*\\.Rdata$"),
        full.names = TRUE
    )
    
    temp <- ls()
    load(rse_path, verbose = verbose)
    obj_name <- ls()[!(ls() %in% c("temp", temp))]
    
    return(get(obj_name))
}

#   Return the colnames in colData(rse) that are not present in that for
#   'rse_ref'
get_unique_names <- function(rse, rse_ref) {
    uniq_names <- colnames(colData(rse))[
        !(colnames(colData(rse)) %in% colnames(colData(rse_ref)))
    ]
    
    return(uniq_names)
}

################################################################################
#   Make sure colData includes the same columns for both objs (FastQC column
#   names are actually experiment-dependent, for example). Adjust year 2's only
################################################################################

#   Load 'rse-gene' objects for all years
rse_gene_1and3 <- load_obj(rse_dir_1and3, 'rse_gene')
rse_gene_2 <- load_obj(rse_dir_2, 'rse_gene')

#   Manually look at the colnames unique to each year
print(get_unique_names(rse_gene_1and3, rse_gene_2))
print(get_unique_names(rse_gene_2, rse_gene_1and3))

#   We see that the columns of year 1 + year 3 are a strict superset of those
#   for year 2. We'll just add NAs for year 2 for those columns
for (missing_name in get_unique_names(rse_gene_1and3, rse_gene_2)) {
    rse_gene_2[[missing_name]] <- NA
}

#   It looks the colData now has the same columns but not in the same order.
#   We'll fix this
colData(rse_gene_2) = colData(rse_gene_2)[
    , match(colnames(colData(rse_gene_1and3)), colnames(colData(rse_gene_2)))
]
stopifnot(
    identical(colnames(colData(rse_gene_1and3)), colnames(colData(rse_gene_2)))
)

################################################################################
#   Subset year 2 objects to use only "main" features (as in year 1, 3) rather
#   than "primary" features. Fix EntrezID
################################################################################

#-------------------------------------------------------------------------------
#   Genes
#-------------------------------------------------------------------------------

#   The "primary" GTF, used for year 2, should contain a strict superset of the
#   genes in the "main" GTF, used for the other years
stopifnot(all(rownames(rse_gene_1and3) %in% rownames(rse_gene_2)))
stopifnot(!all(rownames(rse_gene_2) %in% rownames(rse_gene_1and3)))

#   Subset the year 2 object to contain only the "main" genes seen in year 1, 3
rse_gene_2 <- rse_gene_2[rownames(rse_gene_2) %in% rownames(rse_gene_1and3), ]

stopifnot(identical(ranges(rse_gene_1and3), ranges(rse_gene_2)))

#   EntrezIDs are not identical for the same genes, due to the method for
#   retrieving EntrezID changing. We'll just use year 1 + 3's EntrezIDs
stopifnot(
    !identical(
        rowRanges(rse_gene_1and3)$EntrezID, rowRanges(rse_gene_2)$EntrezID
    )
)
rowRanges(rse_gene_2)$EntrezID <- rowRanges(rse_gene_1and3)$EntrezID

#   Save the updated year 2 object
save(
    rse_gene_2,
    file = file.path(
        out_dir, paste0("rse_gene_n", ncol(rse_gene_2), ".Rdata")
    )
)

metrics = colData(rse_gene_2)

rm(rse_gene_1and3, rse_gene_2)

#-------------------------------------------------------------------------------
#   Exons
#-------------------------------------------------------------------------------

#   Load 'rse-exon' objects for all years
rse_exon_1and3 <- load_obj(rse_dir_1and3, 'rse_exon')
rse_exon_2 <- load_obj(rse_dir_2, 'rse_exon')

#   Ranges of year 2 should be a strict superset of ranges of year 1, 3
stopifnot(
    length(
        ranges(rse_exon_1and3)[!(ranges(rse_exon_1and3) %in% ranges(rse_exon_2))]
    ) == 0
)

stopifnot(
    length(
        ranges(rse_exon_2)[!(ranges(rse_exon_2) %in% ranges(rse_exon_1and3))]
    ) > 0
)

#   Subset year 2 object to have just ranges seen in year 1, 3
rse_exon_2 <- rse_exon_2[ranges(rse_exon_2) %in% ranges(rse_exon_1and3), ]
stopifnot(identical(ranges(rse_exon_1and3), ranges(rse_exon_2)))

#   As we did for the genes, we'll use year 1 EntrezIDs (see above)
rowRanges(rse_exon_2)$EntrezID <- rowRanges(rse_exon_1and3)$EntrezID

#   Use the updated colData
colData(rse_exon_2) = metrics

#   Save the updated year 2 object
save(
    rse_exon_2,
    file = file.path(
        out_dir, paste0("rse_exon_n", ncol(rse_exon_2), ".Rdata")
    )
)

rm(rse_exon_1and3, rse_exon_2)

#-------------------------------------------------------------------------------
#   Transcripts
#-------------------------------------------------------------------------------

#   Load 'rse-tx' objects for all years
rse_tx_1and3 <- load_obj(rse_dir_1and3, 'rse_tx')
rse_tx_2 <- load_obj(rse_dir_2, 'rse_tx')

#   Transcripts are identical
stopifnot(identical(ranges(rse_tx_1and3), ranges(rse_tx_2)))

#   Use the updated colData
colData(rse_tx_2) = metrics

#   Save the updated year 2 object
save(
    rse_tx_2,
    file = file.path(
        out_dir, paste0("rse_tx_n", ncol(rse_tx_2), ".Rdata")
    )
)

rm(rse_tx_1and3, rse_tx_2)
gc()

#-------------------------------------------------------------------------------
#   Junctions
#-------------------------------------------------------------------------------

#   Load 'rse-jx' objects for all years
rse_jx_1and3 <- load_obj(rse_dir_1and3, 'rse_jx')
rse_jx_2 <- load_obj(rse_dir_2, 'rse_jx')

#   We only used the canonical chromosomes for year 1, 3's junctions, but
#   included additional scaffolds for year 2. Subset year 2 to use the same
#   sequences
rse_jx_2 <- rse_jx_2[seqnames(rse_jx_2) %in% seqnames(rse_jx_1and3), ]

#   Get an idea of how ranges overlap between years
length(ranges(rse_jx_1and3)[ranges(rse_jx_1and3) %in% ranges(rse_jx_2)])
length(ranges(rse_jx_2)[ranges(rse_jx_2) %in% ranges(rse_jx_1and3)])

#   Use the updated colData
colData(rse_jx_2) = metrics

#   Save the updated year 2 object
save(
    rse_jx_2,
    file = file.path(
        out_dir, paste0("rse_jx_n", ncol(rse_jx_2), ".Rdata")
    )
)

#   Update the CSV file of metrics
write.csv(
    data.frame(metrics),
    file = file.path(out_dir, "read_and_alignment_metrics.csv")
)

session_info()
