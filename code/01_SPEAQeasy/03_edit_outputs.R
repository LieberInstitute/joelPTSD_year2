#   The colData/ metrics need to be adjusted in year 2 for two different
#   reasons:
#
#   - many columns represent the same information between Y1 and Y2, but have
#     different colnames
#   - the Y2 colData/ metrics don't yet have phenotype information attached, as
#     the Y1 objects do

library('SummarizedExperiment')
library('here')
library('readxl')

rse_path_y1 = '/dcl02/lieber/ajaffe/kleinman_PTSD/count_data/rse_gene_Kleinman_R01_PTSD_n225_Annotated.Rdata'
rse_path_y2 = here(
    'processed-data', '01_SPEAQeasy', 'pipeline_output', 'count_objects',
    'rse_gene_Jlab_experiment_n238.Rdata'
)

#   Not the right path. Does one exist?
y2_meta_path = "/dcl02/lieber/ajaffe/kleinman_PTSD/Joel’s PTSD R01 Y2 RNA Plate Map.xlsx"

################################################################################
#   Functions
################################################################################

#   Rename a column in the colData of 'rse' from 'old_name' to 'new_name'
rename_column = function(rse, old_name, new_name) {
    stopifnot(old_name %in% colnames(colData(rse)))
    stopifnot(! (new_name %in% colnames(colData(rse))))
    
    rse[[new_name]] = rse[[old_name]]
    rse[[old_name]] = NULL
    return(rse)
}

#   Return the colnames in colData(rse) that are not present in that for
#   'rse_ref'
get_unique_names = function(rse, rse_ref) {
    uniq_names = colnames(colData(rse))[
            ! (colnames(colData(rse)) %in% colnames(colData(rse_ref)))
    ]
    
    return(uniq_names)
}

################################################################################
#   First adjust Y2 colData to have compatible colnames with Y1
################################################################################

#   Load 'rse-gene' objects for both years
load(rse_path_y1, verbose = TRUE)
rse_gene_y1 = rse_gene

load(rse_path_y2, verbose = TRUE)
rse_gene_y2 = rse_gene
rm(rse_gene)

#   Manually look at the colnames unique to each year
print(get_unique_names(rse_gene_y1, rse_gene_y2))
print(get_unique_names(rse_gene_y2, rse_gene_y1))

#-------------------------------------------------------------------------------
#   One of the major sources of difference is the naming of FastQC-related
#   columns. Let's adjust these for the Y2 objects to match Y1
#-------------------------------------------------------------------------------

old_names = c(
    "basic_statistics", "per_base_sequence_quality",
    "per_tile_sequence_quality", "per_sequence_quality_scores",
    "per_base_sequence_content", "per_sequence_gc_content",
    "per_base_n_content", "sequence_length_distribution",
    "sequence_duplication_levels", "overrepresented_sequences",
    "adapter_content", "kmer_content"
)

new_names = c(
    "FQCbasicStats", "perBaseQual", "perTileQual", "perSeqQual",
    "perBaseContent", "GCcontent", "Ncontent", "SeqLengthDist",
    "SeqDuplication", "OverrepSeqs", "AdapterContent","KmerContent"
)

stopifnot(length(old_names) == length(new_names))
for (i in 1:length(old_names)) {
    rse_gene_y2 = rename_column(rse_gene_y2, old_names[i], new_names[i])
}

#-------------------------------------------------------------------------------
#   Add columns to each object that are missing from the other
#-------------------------------------------------------------------------------

#   In Y1 and must be added to Y2
missing_cols = c(
    "Adapter88-89_R1", "Adapter88-89_R2", "ERCCsumLogErr",
    "gene_Unassigned_Nonjunction"
)
for (missing_name in missing_cols) {
    rse_gene_y2[[missing_name]] = NA
}

#   In Y2 and must be added to Y1
for (missing_name in get_unique_names(rse_gene_y2, rse_gene_y1)) {
    rse_gene_y1[[missing_name]] = NA
}

#   Look at the columns unique to Y1
print(get_unique_names(rse_gene_y1, rse_gene_y2))

################################################################################
#   Add phenotype data to Y2 object
################################################################################

#   Doesn't currently work
pd_y2 = as.data.frame(read_excel(y2_meta_path))
