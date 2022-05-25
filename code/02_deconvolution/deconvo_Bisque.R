
library("SingleCellExperiment")
library("jaffelab")
library("xbioc")
library("BisqueRNA")
library("dplyr")
library("tidyr")
library("here")
library("sessioninfo")

#### Load Data ####
## Load rse_gene data
load(here("processed-data", "03_SPEAQeasy_Y2", "updated_output", "merged", "rse_gene_pd_n688.Rdata"), verbose = TRUE)
dim(rse_gene)
# [1] 58037   688

## use ensemblID as rownames to match sn Data
rownames(rse_gene) <- rowData(rse_gene)$ensemblID

## sce Data
load("/dcl01/lieber/ajaffe/lab/deconvolution_bsp2/data/sce_pan.v2.Rdata", verbose = TRUE)
load("/dcl01/lieber/ajaffe/lab/deconvolution_bsp2/data/marker_stats_pan.v2.Rdata", verbose = TRUE)

marker_stats %>%
    filter(rank_ratio <= 25, gene %in% rownames(rse_gene)) %>%
    count(cellType.target)

## all genes present
# # A tibble: 10 × 2
#    cellType.target     n
#    <fct>           <int>
#  1 Astro              25
#  2 Endo               25
#  3 Macro              25
#  4 Micro              25
#  5 Mural              25
#  6 Oligo              25
#  7 OPC                25
#  8 Tcell              25
#  9 Excit              25
# 10 Inhib              25

marker_stats2 <- marker_stats %>%
    filter(gene %in% rownames(rse_gene)) %>%
    filter(cellType.target != "Macro" | (cellType.target == "Macro" & rank_ratio <= 23)) %>%
    group_by(cellType.target) %>%
    mutate(rank_ratio = row_number())

marker_stats2 %>%
    filter(rank_ratio <= 25) %>%
    dplyr::count(cellType.target)
# # A tibble: 10 × 2
# # Groups:   cellType.target [10]
#    cellType.target     n
#    <fct>           <int>
#  1 Astro              25
#  2 Endo               25
#  3 Macro              23
#  4 Micro              25
#  5 Mural              25
#  6 Oligo              25
#  7 OPC                25
#  8 Tcell              25
#  9 Excit              25
# 10 Inhib              25

marker_genes <- marker_stats2 %>%
    filter(rank_ratio <= 25) %>%
    pull(gene)

length(marker_genes)
# [1] 248

## save table of marker genes
## marker stats
rd <- as.data.frame(rowData(rse_gene)[marker_genes, ]) %>%
    dplyr::select(ensemblID, gencodeID)

marker_table <- marker_stats2 %>%
    filter(rank_ratio <= 25) %>%
    dplyr::select(ensemblID = gene, `Cell Type` = cellType.target, Symbol) %>%
    left_join(rd) %>%
    dplyr::select(gencodeID, `Cell Type`, Symbol)

readr::write_csv(marker_table, file = here("processed-data", "02_deconvolution", "PTSD_deconvolution_markers.csv"))

#### create expression set ####
exp_set_bulk <- ExpressionSet(
    assayData = assays(rse_gene)$counts[marker_genes, ],
    phenoData = AnnotatedDataFrame(
        as.data.frame(colData(rse_gene))[c("SAMPLE_ID")]
    )
)

exp_set_sce <- ExpressionSet(
    assayData = as.matrix(assays(sce_pan)$counts),
    phenoData = AnnotatedDataFrame(
        as.data.frame(colData(sce_pan))[c("cellType.Broad", "uniqueID", "donor")]
    )
)

## Check for nuclei w/ zero expression in marker genes
exp_set_sce <- exp_set_sce[marker_genes, ]
zero_cell_filter <- colSums(exprs(exp_set_sce)) != 0
message("Exclude ", sum(!zero_cell_filter), " cells")
# Exclude 2 cells

exp_set_sce <- exp_set_sce[, zero_cell_filter]

#### Run Bisque ####
est_prop <- ReferenceBasedDecomposition(
    bulk.eset = exp_set_bulk,
    sc.eset = exp_set_sce,
    cell.types = "cellType.Broad",
    subject.names = "donor",
    use.overlap = FALSE
)

est_prop$bulk.props <- t(est_prop$bulk.props)

## Save proportions as csv
write.csv(est_prop$bulk.props, file = here("processed-data", "02_deconvolution", "PTSD_deconvolution_results.csv"))

est_prop$Est.prop.long <- est_prop$bulk.props %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Sample") %>%
    pivot_longer(!Sample, names_to = "cell_type", values_to = "prop")

## Add long data and save
round(colMeans(est_prop$bulk.props), 3)
# Astro  Endo Macro Micro Mural Oligo   OPC Tcell Excit Inhib
# 0.066 0.003 0.005 0.047 0.010 0.351 0.062 0.008 0.094 0.353

summary(est_prop$bulk.props)
#     Astro              Endo               Macro              Micro
# Min.   :0.00000   Min.   :0.0000000   Min.   :0.000000   Min.   :0.00000
# 1st Qu.:0.04701   1st Qu.:0.0000000   1st Qu.:0.000000   1st Qu.:0.03268
# Median :0.06590   Median :0.0008699   Median :0.001698   Median :0.04870
# Mean   :0.06634   Mean   :0.0031168   Mean   :0.005302   Mean   :0.04663
# 3rd Qu.:0.08629   3rd Qu.:0.0059085   3rd Qu.:0.009904   3rd Qu.:0.06286
# Max.   :0.18017   Max.   :0.0214793   Max.   :0.038200   Max.   :0.11878
#     Mural              Oligo             OPC              Tcell
# Min.   :0.000000   Min.   :0.0000   Min.   :0.00000   Min.   :0.000000
# 1st Qu.:0.000000   1st Qu.:0.2320   1st Qu.:0.03806   1st Qu.:0.000000
# Median :0.003759   Median :0.3265   Median :0.05703   Median :0.003830
# Mean   :0.009970   Mean   :0.3511   Mean   :0.06235   Mean   :0.008261
# 3rd Qu.:0.017108   3rd Qu.:0.4399   3rd Qu.:0.08048   3rd Qu.:0.014622
# Max.   :0.064820   Max.   :1.0000   Max.   :0.19508   Max.   :0.045075
#     Excit             Inhib
# Min.   :0.00000   Min.   :0.0000
# 1st Qu.:0.05674   1st Qu.:0.2808
# Median :0.08758   Median :0.3429
# Mean   :0.09429   Mean   :0.3526
# 3rd Qu.:0.13797   3rd Qu.:0.4227
# Max.   :0.23441   Max.   :1.0000

## Explore a bit more the samples that have no neurons
with(as.data.frame(est_prop$bulk.props), addmargins(table(
    "No excitatory" = Excit == 0, "No inhibitory" = Inhib == 0
)))
#              No inhibitory
# No excitatory FALSE TRUE Sum
#         FALSE   422    0 422
#         TRUE     39    2  41
#         Sum     461    2 463

stopifnot(identical(rownames(est_prop$bulk.props), rse_gene$RNum))

## Explore the samples with no excitatory neurons
rse_gene$no_excit <- est_prop$bulk.props[, "Excit"] == 0
no_excit <- as.data.frame(subset(colData(rse_gene)[, c(
    "Age",
    "Region",
    "Sex",
    "Dx",
    "Race",
    "Dataset",
    "RIN",
    "no_excit"
)], no_excit == TRUE))
summary(data.frame(lapply(no_excit, function(x) {
    if (is.character(x)) {
        return(factor(x))
    } else {
        return(x)
    }
})))
#      Age                  Region   Sex          Dx       Race
# Min.   :17.47   Central Amyg :19   F:20   Control:20   AA  : 4
# 1st Qu.:31.38   Dentate Gyrus:17   M:21   MDD    :10   AS  : 1
# Median :45.13   mPFC         : 5          PTSD   :11   CAUC:35
# Mean   :43.56                                          HISP: 1
# 3rd Qu.:54.71
# Max.   :69.03
#       Dataset        RIN        no_excit
# PTSD_Year1:18   Min.   :5.100   Mode:logical
# PTSD_Year2:23   1st Qu.:6.000   TRUE:41
#                 Median :6.300
#                 Mean   :6.395
#                 3rd Qu.:6.600
#                 Max.   :9.700

## Hm... I don't see any obvious patterns. Maybe the brain region?

## Save full results
est_prop_bisque <- est_prop
save(est_prop_bisque, file = here("processed-data", "02_deconvolution", "est_prop_Bisque.Rdata"))

# sgejobs::job_single('deconvo_Bisque', create_shell = TRUE, queue= 'bluejay', memory = '40G', command = "Rscript deconvo_Bisque.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
