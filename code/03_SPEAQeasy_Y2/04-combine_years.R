#   For convenience, we'd like to combine RangedSummarizedExperiment objects
#   across years 1-3. Objects were adjusted as needed for compatibility in the
#   last script, 03-edit_outputs.R, to make combining straightforward

library("SummarizedExperiment")
library("here")
library("sessioninfo")

rse_dir_1and3 <- here(
    "processed-data", "01_SPEAQeasy", "updated_output", "merged"
)
rse_dir_2 <- here(
    "processed-data", "03_SPEAQeasy_Y2", "updated_output", "year_2"
)

out_dir <- here("processed-data", "03_SPEAQeasy_Y2", "updated_output", "merged")
dir.create(out_dir, showWarnings = FALSE)

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

################################################################################
#   Combine genes, exons, and transcripts
################################################################################

for (var_name in c('rse_gene', 'rse_exon', 'rse_tx')) {
    obj_1and3 <- load_obj(rse_dir_1and3, var_name)
    obj_2 <- load_obj(rse_dir_2, var_name)
    
    if (var_name %in% c('rse_gene', 'rse_exon')) {
        #   Recompute meanExprs by taking into account all samples across the 3
        #   years
        temp <- (rowRanges(obj_1and3)$meanExprs * ncol(obj_1and3) +
                     rowRanges(obj_2)$meanExprs * ncol(obj_2)) /
            (ncol(obj_1and3) + ncol(obj_2))
        rowRanges(obj_1and3)$meanExprs <- temp
        rowRanges(obj_2)$meanExprs <- temp
    }
    
    #   No differences should now exists between rowRanges
    stopifnot(identical(rowRanges(obj_1and3), rowRanges(obj_1and3)))
    
    #   Combine and save
    obj_combined = cbind(obj_1and3, obj_2)

    assign(var_name, obj_combined)
    save(
        list = var_name,
        file = file.path(
            out_dir,
            paste0(var_name, "_n", ncol(obj_combined), ".Rdata")
        )
    )
    
    #   Explicitly free memory
    rm(obj_1and3, obj_2, obj_combined)
    rm(list= var_name)
    gc()
}

################################################################################
#   Combine junctions
################################################################################

#   Load 'rse-jx' objects for all years
rse_jx_1and3 <- load_obj(rse_dir_1and3, 'rse_jx')
rse_jx_2 <- load_obj(rse_dir_2, 'rse_jx')

#-------------------------------------------------------------------------------
#   rowRanges
#-------------------------------------------------------------------------------

#   Re-compute meanExprs column for junctions seen in both year 1,3 and year 2
shared_names <- names(ranges(rse_jx_1and3))[
    names(ranges(rse_jx_1and3)) %in% names(ranges(rse_jx_2))
]

mean_exprs_1and3 <- rowRanges(rse_jx_1and3)$meanExprs[shared_names]
mean_exprs_2 <- rowRanges(rse_jx_2)$meanExprs[shared_names]

shared_mean <- (mean_exprs_1and3 * ncol(rse_jx_1and3) +
    mean_exprs_2 * ncol(rse_jx_2)) /
    (ncol(rse_jx_1and3) + ncol(rse_jx_2))

#   Take the union of rowRanges
gr_union <- unique(c(rowRanges(rse_jx_1and3), rowRanges(rse_jx_2)))
gr_union$meanExprs[shared_names] <- shared_mean

#-------------------------------------------------------------------------------
#   assays
#-------------------------------------------------------------------------------

#   Initialize NA counts assays for the future merged rse_jx
counts <- matrix(
    NA,
    nrow = length(gr_union), ncol = ncol(rse_jx_1and3) + ncol(rse_jx_2),
    dimnames = list(
        names(gr_union), c(colnames(rse_jx_1and3), colnames(rse_jx_2))
    )
)

#   Fill counts with the assays from each year
counts[rownames(rse_jx_1and3), colnames(rse_jx_1and3)] <- assays(rse_jx_1and3)$counts
counts[rownames(rse_jx_2), colnames(rse_jx_2)] <- assays(rse_jx_2)$counts

#-------------------------------------------------------------------------------
#   Construct combined 'rse_jx' SummarizedExperiment
#-------------------------------------------------------------------------------

rse_jx <- SummarizedExperiment(
    assays = list("counts" = counts), rowRanges = gr_union,
    colData = rbind(colData(rse_jx_1and3), colData(rse_jx_2))
)

save(
    rse_jx,
    file = here(out_dir, "merged", paste0("rse_jx_n", ncol(rse_jx), ".Rdata"))
)

session_info()
