
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
load(here("processed-data","01_SPEAQeasy","updated_output", "merged", "rse_gene.y1y2_n463.Rdata"), verbose = TRUE)
## use ensemblID as rownames to match sn Data
rownames(rse_gene) <- rowData(rse_gene)$ensemblID

## sce Data
load("/dcl01/lieber/ajaffe/lab/deconvolution_bsp2/data/sce_pan.v2.Rdata", verbose = TRUE)
load("/dcl01/lieber/ajaffe/lab/deconvolution_bsp2/data/marker_stats_pan.v2.Rdata", verbose = TRUE)

marker_stats %>%
  filter(rank_ratio <= 25, gene %in% rownames(rse_gene)) %>%
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
  filter(gene %in% rownames(rse_gene)) %>%
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
rd <- as.data.frame(rowData(rse_gene)[marker_genes,]) %>%
  dplyr::select(ensemblID, gencodeID)

marker_table <- marker_stats2 %>% filter(rank_ratio <= 25) %>%
  dplyr::select(ensemblID = gene, `Cell Type` = cellType.target, Symbol) %>%
  left_join(rd) %>%
  dplyr::select(gencodeID, `Cell Type`, Symbol)

write_csv(marker_table, file = here("processed-data", "02_deconvolution","PTSD_deconvolution_markers.csv"))

#### create expression set ####
exp_set_bulk <- ExpressionSet(assayData = assays(rse_gene)$counts[marker_genes,],
                              phenoData=AnnotatedDataFrame(
                                as.data.frame(colData(rse_gene))[c("SAMPLE_ID")]))

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

# sgejobs::job_single('deconvo_Bisque', create_shell = TRUE, queue= 'bluejay', memory = '20G', command = "Rscript deconvo_Bisque.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
