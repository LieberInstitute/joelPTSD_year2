#   For convenience, we'd like to combine RangedSummarizedExperiment objects
#   across year 1 and year 2. There are subtle differences due to code changes
#   made between processing dates, which are resolved here before combining the
#   objects.

library('SummarizedExperiment')
library('here')
library('sessioninfo')

rse_dir = here('processed-data', '01_SPEAQeasy', 'updated_output')

load_obj = function(var_name, year, verbose = TRUE) {
    this_dir = here(rse_dir, paste0('year_', year))
    rse_path = list.files(
        this_dir,
        pattern = paste0('^', var_name, '_n.*\\.Rdata$'),
        full.names = TRUE
    )
    
    temp = ls()
    load(rse_path, verbose = verbose)
    obj_name = ls()[!(ls() %in% c('temp', temp))]
    
    return(get(obj_name))
}

################################################################################
#   Genes
################################################################################

rse_gene_y1 = load_obj('rse_gene', 1)
rse_gene_y2 = load_obj('rse_gene', 2)

#   The "primary" GTF, used for year 2, should contain a strict superset of the
#   genes in the "main" GTF, used for year 1.
stopifnot(all(rownames(rse_gene_y1) %in% rownames(rse_gene_y2)))
stopifnot(!all(rownames(rse_gene_y2) %in% rownames(rse_gene_y1)))

#   Subset the year 2 object to contain only the "main" genes seen in year 1
rse_gene_y2 = rse_gene_y2[rownames(rse_gene_y2) %in% rownames(rse_gene_y1),]

stopifnot(all(rownames(rse_gene_y1) == rownames(rse_gene_y2)))

#   Assess why ranges for the same genes differ between years
equal_genes = ranges(rse_gene_y1) == ranges(rse_gene_y2)
length(which(equal_genes)) / length(equal_genes)
which(!equal_genes)

#   Ranges only differ by end coordinate...
stopifnot(all(start(ranges(rse_gene_y1)) == start(ranges(rse_gene_y2))))

#   Looking up a few genes at random, it looks like the year 1 coordinates
#   are always accurate, but some year 2 end coordinates are wrong. For now,
#   just use the year 1 coords
ranges(rse_gene_y2) = ranges(rse_gene_y1)

#   Now it looks like 'EntrezID' columns contain different information for some
#   rows.
equal_entrez = rowRanges(rse_gene_y1)$EntrezID == rowRanges(rse_gene_y2)$EntrezID
equal_entrez = is.na(equal_entrez) | equal_entrez

#   Somewhat surprisingly, a difference in EntrezID between years is not always
#   due to differences in the ranges. The difference is in some cases likely
#   due to the method for retrieving EntrezID changing. The method matters
#   because more than 1 EntrezID can exist for a single gene
table(which(!equal_entrez) %in% which(!equal_genes))

#   We'll just use year 1's EntrezIDs, especially since year 2 has some bad
#   ranges that might result in bad EntrezIDs
rowRanges(rse_gene_y2)$EntrezID = rowRanges(rse_gene_y1)$EntrezID

#   Save the adjusted copy of 'rse_gene_y2'
dir.create(here(rse_dir, 'year_2', 'adjusted_ranges'), showWarnings = FALSE)
out_file = here(
    rse_dir, 'year_2', 'adjusted_ranges', 
    paste0('rse_gene_n', ncol(rse_gene_y2), '.Rdata')
)
save(rse_gene_y2, file = out_file)

#   Now the 'meanExprs' column is different (which should actually be expected);
#   recompute meanExprs by taking into account all samples
temp = (rowRanges(rse_gene_y1)$meanExprs * ncol(rse_gene_y1) + 
    rowRanges(rse_gene_y2)$meanExprs * ncol(rse_gene_y2)) / 
    (ncol(rse_gene_y1) + ncol(rse_gene_y2))
rowRanges(rse_gene_y1)$meanExprs = temp
rowRanges(rse_gene_y2)$meanExprs = temp

#   Finally we can 'cbind' the objects as planned
rse_gene = cbind(rse_gene_y1, rse_gene_y2)

#   Save the merged 'rse_gene'
dir.create(here(rse_dir, 'merged'), showWarnings = FALSE)
out_file = here(
    rse_dir, 'merged', paste0('rse_gene_n', ncol(rse_gene), '.Rdata')
)
save(rse_gene, file = out_file)

################################################################################
#   Exons
################################################################################

rse_exon_y1 = load_obj('rse_exon', 1)
rse_exon_y2 = load_obj('rse_exon', 2)

#   Ranges of year 2 should be a strict superset of ranges of year 1
stopifnot(
    length(
        ranges(rse_exon_y1)[!(ranges(rse_exon_y1) %in% ranges(rse_exon_y2))]
    ) == 0
)

stopifnot(
    length(
        ranges(rse_exon_y2)[!(ranges(rse_exon_y2) %in% ranges(rse_exon_y1))]
    ) > 0
)

#   Subset year 2 object to have just ranges seen in year 1
rse_exon_y2 = rse_exon_y2[ranges(rse_exon_y2) %in% ranges(rse_exon_y1),]

#   As we did for the genes, we'll use year 1 EntrezIDs (see above)
rowRanges(rse_exon_y2)$EntrezID = rowRanges(rse_exon_y1)$EntrezID

#   Save the adjusted copy of 'rse_exon_y2'
out_file = here(
    rse_dir, 'year_2', 'adjusted_ranges', 
    paste0('rse_exon_n', ncol(rse_exon_y2), '.Rdata')
)
save(rse_exon_y2, file = out_file)

#   Again, recompute 'meanExpr' since we're combining data
#   Now the 'meanExprs' column is different (which should actually be expected);
#   recompute meanExprs by taking into account all samples
temp = (rowRanges(rse_exon_y1)$meanExprs * ncol(rse_exon_y1) + 
            rowRanges(rse_exon_y2)$meanExprs * ncol(rse_exon_y2)) / 
    (ncol(rse_exon_y1) + ncol(rse_exon_y2))
rowRanges(rse_exon_y1)$meanExprs = temp
rowRanges(rse_exon_y2)$meanExprs = temp

#   We can now 'cbind' the objects as planned
rse_exon = cbind(rse_exon_y1, rse_exon_y2)

#   Save the merged 'rse_exon'
out_file = here(
    rse_dir, 'merged', paste0('rse_exon_n', ncol(rse_exon), '.Rdata')
)
save(rse_exon, file = out_file)

################################################################################
#   Transcripts
################################################################################

rse_tx_y1 = load_obj('rse_tx', 1)
rse_tx_y2 = load_obj('rse_tx', 2)

#   Transcripts are identical
stopifnot(identical(ranges(rse_tx_y1), ranges(rse_tx_y2)))

#   The 'counts' assay is only defined for year 2. Just make it NA in year 1
assays(rse_tx_y1, withDimnames=FALSE)$counts = matrix(
    rep(NA, dim(rse_tx_y1)[1] * dim(rse_tx_y1)[2]),
    nrow = dim(rse_tx_y1)[1]
)

#   We can now 'cbind' the objects as planned
rse_tx = cbind(rse_tx_y1, rse_tx_y2)

#   Save the merged 'rse_tx'
out_file = here(
    rse_dir, 'merged', paste0('rse_tx_n', ncol(rse_tx), '.Rdata')
)
save(rse_tx, file = out_file)

################################################################################
#   Junctions
################################################################################

rse_jx_y1 = load_obj('rse_jx', 1)
rse_jx_y2 = load_obj('rse_jx', 2)

#   We only used the canonical chromosomes for year 1's junctions, but included
#   additional scaffolds for year 2. Subset year 2 to use the same sequences
rse_jx_y2 = rse_jx_y2[seqnames(rse_jx_y2) %in% seqnames(rse_jx_y1),]

#   Get an idea of how ranges overlap between years
length(ranges(rse_jx_y1)[ranges(rse_jx_y1) %in% ranges(rse_jx_y2)])
length(ranges(rse_jx_y2)[ranges(rse_jx_y2) %in% ranges(rse_jx_y1)])

#-------------------------------------------------------------------------------
#   rowRanges
#-------------------------------------------------------------------------------

#   Re-compute meanExprs column for junctions seen in both year 1 and year 2
shared_names = names(ranges(rse_jx_y1))[
    names(ranges(rse_jx_y1)) %in% names(ranges(rse_jx_y2))
]

mean_exprs_y1 = rowRanges(rse_jx_y1)$meanExprs[shared_names]
mean_exprs_y2 = rowRanges(rse_jx_y2)$meanExprs[shared_names]

shared_mean = (mean_exprs_y1 * ncol(rse_jx_y1) +
    mean_exprs_y2 * ncol(rse_jx_y2)) / 
    (ncol(rse_jx_y1) + ncol(rse_jx_y2))

#   Take the union of rowRanges
gr_union = unique(c(rowRanges(rse_jx_y1), rowRanges(rse_jx_y2)))
gr_union$meanExprs[shared_names] = shared_mean

#-------------------------------------------------------------------------------
#   assays
#-------------------------------------------------------------------------------

#   Initialize NA counts assays for the future merged rse_jx
counts = matrix(
    NA, nrow = length(gr_union), ncol = ncol(rse_jx_y1) + ncol(rse_jx_y2),
    dimnames = list(
        names(gr_union), c(colnames(rse_jx_y1), colnames(rse_jx_y2))
    )
)

#   Fill counts with the assays from each year
counts[rownames(rse_jx_y1), colnames(rse_jx_y1)] = assays(rse_jx_y1)$counts
counts[rownames(rse_jx_y2), colnames(rse_jx_y2)] = assays(rse_jx_y2)$counts

#-------------------------------------------------------------------------------
#   Construct combined 'rse_jx' SummarizedExperiment
#-------------------------------------------------------------------------------

rse_jx = SummarizedExperiment(
    assays = list('counts' = counts), rowRanges = gr_union,
    colData = rbind(colData(rse_jx_y1), colData(rse_jx_y2))
)

#   Save
out_file = here(
    rse_dir, 'merged', paste0('rse_jx_n', ncol(rse_jx), '.Rdata')
)
save(rse_jx, file = out_file)

session_info()
