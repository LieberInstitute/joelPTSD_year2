
library("SummarizedExperiment")
library("SingleCellExperiment")
library("jaffelab")
library("xbioc")
library("BisqueRNA")
library("tidyverse")
library("here")
library("sessioninfo")

#### Load Data ####
## Load rse_gene data
load(here("processed-data","01_SPEAQeasy","updated_output", "year_2", "adjusted_ranges", "rse_gene_n238.Rdata"), verbose = TRUE)
## use ensemblID as rownames to match sn Data
rownames(rse_gene_y2) <- rowData(rse_gene_y2)$ensemblID

## sce Data
load("/dcl01/lieber/ajaffe/lab/deconvolution_bsp2/data/sce_pan.v2.Rdata", verbose = TRUE)
load("/dcl01/lieber/ajaffe/lab/deconvolution_bsp2/data/marker_stats_pan.v2.Rdata", verbose = TRUE)

marker_stats %>%
  filter(rank_ratio <= 25, gene %in% rownames(rse_gene_y2)) %>%
  count(cellType.target)

## all genes present
# cellType.target     n
# <fct>           <int>
# 1 Astro              25
# 2 Endo               25
# 3 Macro              25
# 4 Micro              25
# 5 Mural              25
# 6 Oligo              25
# 7 OPC                25
# 8 Tcell              25
# 9 Excit              25
# 10 Inhib              25

marker_stats2 <- marker_stats %>%
  filter(gene %in% rownames(rse_gene_y2)) %>%
  filter(cellType.target != "Macro" | (cellType.target == "Macro" & rank_ratio <= 23)) %>%
  group_by(cellType.target) %>%
  mutate(rank_ratio = row_number())

marker_stats2%>%
  filter(rank_ratio <= 25) %>% 
  dplyr::count(cellType.target)

marker_genes <- marker_stats2 %>%
  filter(rank_ratio <= 25) %>%
  pull(gene)

length(marker_genes)
# [1] 248

## save table of marker genes
## marker stats
rd <- as.data.frame(rowData(rse_gene_y2)[marker_genes,]) %>% 
  dplyr::select(ensemblID, gencodeID)

marker_table <- marker_stats2 %>% filter(rank_ratio <= 25) %>%
  dplyr::select(ensemblID = gene, `Cell Type` = cellType.target, Symbol) %>%
  left_join(rd) %>%
  dplyr::select(gencodeID, `Cell Type`, Symbol)

write_csv(marker_table, file = here("processed-data", "02_deconvolution","PTSD_deconvolution_markers.csv"))

#### create expression set ####
exp_set_bulk <- ExpressionSet(assayData = assays(rse_gene_y2)$counts[marker_genes,],
                              phenoData=AnnotatedDataFrame(
                                as.data.frame(colData(rse_gene_y2))[c("SAMPLE_ID")]))

exp_set_sce <- ExpressionSet(assayData = as.matrix(assays(sce_pan)$counts),
                             phenoData=AnnotatedDataFrame(
                               as.data.frame(colData(sce_pan))[c("cellType.Broad", "uniqueID","donor")]))

## Check for nuclei w/ zero expression in marker genes
exp_set_sce <- exp_set_sce[marker_genes,]
zero_cell_filter <- colSums(exprs(exp_set_sce)) != 0
message("Exclude ",sum(!zero_cell_filter), " cells")
# Exclude 2 cells

exp_set_sce <- exp_set_sce[,zero_cell_filter]

#### Run Bisque ####
est_prop <- ReferenceBasedDecomposition(bulk.eset = exp_set_bulk,
                                        sc.eset = exp_set_sce,
                                        cell.types = "cellType.Broad",
                                        subject.names = "donor",
                                        use.overlap = FALSE)

est_prop$bulk.props <- t(est_prop$bulk.props)

## Save proportions as csv
write.csv(est_prop$bulk.props, file = here("processed-data", "02_deconvolution","PTSD_deconvolution_results.csv"))

est_prop$Est.prop.long <- est_prop$bulk.props %>%
  as.data.frame() %>%
  rownames_to_column("Sample")%>%
  pivot_longer(!Sample, names_to = "cell_type", values_to = "prop")

## Add long data and save
round(colMeans(est_prop$bulk.props), 3)
# Astro  Endo Macro Micro Mural Oligo   OPC Tcell Excit Inhib 
# 0.066 0.003 0.005 0.047 0.009 0.351 0.063 0.008 0.095 0.353 

summary(est_prop$bulk.props)

## Save full results
est_prop_bisque <- est_prop
save(est_prop_bisque, file = here("processed-data", "02_deconvolution","est_prop_Bisque.Rdata"))

# sgejobs::job_single('deconvo_Bisque', create_shell = TRUE, queue= 'bluejay', memory = '10G', command = "Rscript deconvo_Bisque.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R Under development (unstable) (2021-11-06 r81149)
# os       CentOS Linux 7 (Core)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2022-05-09
# pandoc   2.11.0.4 @ /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-devel/bin/pandoc
# 
# ─ Packages ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version  date (UTC) lib source
# AnnotationDbi        * 1.58.0   2022-04-26 [2] Bioconductor
# assertthat             0.2.1    2019-03-21 [2] CRAN (R 4.1.0)
# backports              1.4.1    2021-12-13 [2] CRAN (R 4.2.0)
# Biobase              * 2.56.0   2022-04-26 [2] Bioconductor
# BiocGenerics         * 0.42.0   2022-04-26 [2] Bioconductor
# BiocManager            1.30.17  2022-04-22 [2] CRAN (R 4.2.0)
# Biostrings             2.64.0   2022-04-26 [2] Bioconductor
# BisqueRNA            * 1.0.5    2021-05-23 [1] CRAN (R 4.2.0)
# bit                    4.0.4    2020-08-04 [2] CRAN (R 4.1.0)
# bit64                  4.0.5    2020-08-30 [2] CRAN (R 4.1.0)
# bitops                 1.0-7    2021-04-24 [2] CRAN (R 4.2.0)
# blob                   1.2.3    2022-04-10 [2] CRAN (R 4.2.0)
# broom                  0.8.0    2022-04-13 [2] CRAN (R 4.2.0)
# cachem                 1.0.6    2021-08-19 [2] CRAN (R 4.2.0)
# callr                  3.7.0    2021-04-20 [2] CRAN (R 4.2.0)
# cellranger             1.1.0    2016-07-27 [2] CRAN (R 4.1.0)
# checkmate              2.1.0    2022-04-21 [2] CRAN (R 4.2.0)
# cli                    3.3.0    2022-04-25 [2] CRAN (R 4.2.0)
# codetools              0.2-18   2020-11-04 [3] CRAN (R 4.2.0)
# colorout             * 1.2-2    2022-04-06 [1] Github (jalvesaq/colorout@79931fd)
# colorspace             2.0-3    2022-02-21 [2] CRAN (R 4.2.0)
# crayon                 1.5.1    2022-03-26 [2] CRAN (R 4.2.0)
# curl                   4.3.2    2021-06-23 [2] CRAN (R 4.2.0)
# DBI                    1.1.2    2021-12-20 [2] CRAN (R 4.2.0)
# dbplyr                 2.1.1    2021-04-06 [2] CRAN (R 4.1.0)
# DelayedArray           0.22.0   2022-04-26 [2] Bioconductor
# digest                 0.6.29   2021-12-01 [2] CRAN (R 4.2.0)
# dplyr                * 1.0.9    2022-04-28 [2] CRAN (R 4.2.0)
# ellipsis               0.3.2    2021-04-29 [2] CRAN (R 4.2.0)
# fansi                  1.0.3    2022-03-24 [2] CRAN (R 4.2.0)
# fastmap                1.1.0    2021-01-25 [2] CRAN (R 4.1.0)
# forcats              * 0.5.1    2021-01-27 [2] CRAN (R 4.1.0)
# fs                     1.5.2    2021-12-08 [2] CRAN (R 4.2.0)
# gargle                 1.2.0    2021-07-02 [2] CRAN (R 4.2.0)
# generics               0.1.2    2022-01-31 [2] CRAN (R 4.2.0)
# GenomeInfoDb         * 1.32.1   2022-04-28 [2] Bioconductor
# GenomeInfoDbData       1.2.8    2022-04-16 [2] Bioconductor
# GenomicRanges        * 1.48.0   2022-04-26 [2] Bioconductor
# ggplot2              * 3.3.6    2022-05-03 [2] CRAN (R 4.2.0)
# glue                   1.6.2    2022-02-24 [2] CRAN (R 4.2.0)
# googledrive            2.0.0    2021-07-08 [2] CRAN (R 4.2.0)
# gtable                 0.3.0    2019-03-25 [2] CRAN (R 4.1.0)
# haven                  2.5.0    2022-04-15 [2] CRAN (R 4.2.0)
# here                 * 1.0.1    2020-12-13 [1] CRAN (R 4.2.0)
# hms                    1.1.1    2021-09-26 [2] CRAN (R 4.2.0)
# httr                   1.4.3    2022-05-04 [2] CRAN (R 4.2.0)
# IRanges              * 2.30.0   2022-04-26 [2] Bioconductor
# jaffelab             * 0.99.32  2022-04-06 [1] Github (LieberInstitute/jaffelab@7b7afe3)
# jsonlite               1.8.0    2022-02-22 [2] CRAN (R 4.2.0)
# KEGGREST               1.36.0   2022-04-26 [2] Bioconductor
# lattice                0.20-45  2021-09-22 [3] CRAN (R 4.2.0)
# lifecycle              1.0.1    2021-09-24 [2] CRAN (R 4.2.0)
# limma                  3.52.0   2022-04-26 [2] Bioconductor
# limSolve               1.5.6    2019-11-14 [1] CRAN (R 4.2.0)
# lpSolve                5.6.15   2020-01-24 [1] CRAN (R 4.2.0)
# lubridate              1.8.0    2021-10-07 [2] CRAN (R 4.2.0)
# magrittr               2.0.3    2022-03-30 [2] CRAN (R 4.2.0)
# MASS                   7.3-57   2022-04-22 [3] CRAN (R 4.2.0)
# Matrix                 1.4-1    2022-03-23 [3] CRAN (R 4.2.0)
# MatrixGenerics       * 1.8.0    2022-04-26 [2] Bioconductor
# matrixStats          * 0.62.0   2022-04-19 [2] CRAN (R 4.2.0)
# memoise                2.0.1    2021-11-26 [2] CRAN (R 4.2.0)
# modelr                 0.1.8    2020-05-19 [2] CRAN (R 4.1.0)
# munsell                0.5.0    2018-06-12 [2] CRAN (R 4.1.0)
# pillar                 1.7.0    2022-02-01 [2] CRAN (R 4.2.0)
# pkgbuild               1.3.1    2021-12-20 [2] CRAN (R 4.2.0)
# pkgconfig              2.0.3    2019-09-22 [2] CRAN (R 4.1.0)
# pkgmaker               0.32.2   2020-10-20 [1] CRAN (R 4.2.0)
# plyr                   1.8.7    2022-03-24 [2] CRAN (R 4.2.0)
# png                    0.1-7    2013-12-03 [2] CRAN (R 4.1.0)
# prettyunits            1.1.1    2020-01-24 [2] CRAN (R 4.1.0)
# processx               3.5.3    2022-03-25 [2] CRAN (R 4.2.0)
# ps                     1.7.0    2022-04-23 [2] CRAN (R 4.2.0)
# purrr                * 0.3.4    2020-04-17 [2] CRAN (R 4.1.0)
# quadprog               1.5-8    2019-11-20 [2] CRAN (R 4.1.0)
# R6                     2.5.1    2021-08-19 [2] CRAN (R 4.2.0)
# rafalib              * 1.0.0    2015-08-09 [1] CRAN (R 4.2.0)
# RColorBrewer           1.1-3    2022-04-03 [2] CRAN (R 4.2.0)
# Rcpp                   1.0.8.3  2022-03-17 [2] CRAN (R 4.2.0)
# RCurl                  1.98-1.6 2022-02-08 [2] CRAN (R 4.2.0)
# readr                * 2.1.2    2022-01-30 [2] CRAN (R 4.2.0)
# readxl                 1.4.0    2022-03-28 [2] CRAN (R 4.2.0)
# registry               0.5-1    2019-03-05 [2] CRAN (R 4.1.0)
# remotes                2.4.2    2021-11-30 [2] CRAN (R 4.2.0)
# reprex                 2.0.1    2021-08-05 [2] CRAN (R 4.2.0)
# reshape2               1.4.4    2020-04-09 [2] CRAN (R 4.1.0)
# rlang                  1.0.2    2022-03-04 [2] CRAN (R 4.2.0)
# rprojroot              2.0.3    2022-04-02 [2] CRAN (R 4.2.0)
# RSQLite                2.2.14   2022-05-07 [2] CRAN (R 4.2.0)
# rstudioapi             0.13     2020-11-12 [2] CRAN (R 4.1.0)
# rvest                  1.0.2    2021-10-16 [2] CRAN (R 4.2.0)
# S4Vectors            * 0.34.0   2022-04-26 [2] Bioconductor
# scales                 1.2.0    2022-04-13 [2] CRAN (R 4.2.0)
# segmented              1.5-0    2022-04-11 [1] CRAN (R 4.2.0)
# sessioninfo          * 1.2.2    2021-12-06 [2] CRAN (R 4.2.0)
# SingleCellExperiment * 1.18.0   2022-04-26 [2] Bioconductor
# stringi                1.7.6    2021-11-29 [2] CRAN (R 4.2.0)
# stringr              * 1.4.0    2019-02-10 [2] CRAN (R 4.1.0)
# SummarizedExperiment * 1.26.1   2022-04-29 [2] Bioconductor
# tibble               * 3.1.7    2022-05-03 [2] CRAN (R 4.2.0)
# tidyr                * 1.2.0    2022-02-01 [2] CRAN (R 4.2.0)
# tidyselect             1.1.2    2022-02-21 [2] CRAN (R 4.2.0)
# tidyverse            * 1.3.1    2021-04-15 [2] CRAN (R 4.2.0)
# tzdb                   0.3.0    2022-03-28 [2] CRAN (R 4.2.0)
# utf8                   1.2.2    2021-07-24 [2] CRAN (R 4.2.0)
# vctrs                  0.4.1    2022-04-13 [2] CRAN (R 4.2.0)
# vroom                  1.5.7    2021-11-30 [2] CRAN (R 4.2.0)
# withr                  2.5.0    2022-03-03 [2] CRAN (R 4.2.0)
# xbioc                * 0.1.19   2022-05-09 [1] Github (renozao/xbioc@1354168)
# xml2                   1.3.3    2021-11-30 [2] CRAN (R 4.2.0)
# xtable                 1.8-4    2019-04-21 [2] CRAN (R 4.1.0)
# XVector                0.36.0   2022-04-26 [2] Bioconductor
# zlibbioc               1.42.0   2022-04-26 [2] Bioconductor
# 
# [1] /users/lhuuki/R/devel
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-devel/R/devel/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-devel/R/devel/lib64/R/library
# 
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
