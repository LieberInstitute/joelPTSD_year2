Deconvolution Analysis
========

Louise A. Huuki-Myers
5/9/2022

### Details
[`deconvo_Bisque.R`](deconvo_Bisque.R): Run deconvolution analysis with Bisque 
([Jew et al., Nature Communications, 2013](https://doi.org/10.1038/s41467-020-15816-6))  

* Used 5 region "pan" brain data from [Tran, Maynard et al., Neuron, 2021](https://github.com/LieberInstitute/10xPilot_snRNAseq-human) 
as snRNA-seq reference  
* ten "broad" cell types  
* top 25 "mean ratio" marker genes (except 23 for Macro), recorded in [`processed-data/02_deconvolution/PTSD_deconvolution_markers.csv`](../../processed-data/02_deconvolution/PTSD_deconvolution_markers.csv)  
* results recorded in [`processed-data/02_deconvolution/PTSD_deconvolution_results.csv`](../../processed-data/02_deconvolution/PTSD_deconvolution_results.csv) & [`processed-data/02_deconvolution/est_prop_Bisque.Rdata`](../../processed-data/02_deconvolution/est_prop_Bisque.Rdata)  

| Cell Type   | Acronym |
| ----------- | ----------- |
| Astrocyte | Astro|
| Endothelial Cells | Endo |
| Macrophage | Macro |
| Microglia | Micro |
| Mural Cells | Mural | 
| Oligodendrocytes | Oligo | 
| Oligodendrocyte Progenitor Cells | OPC | 
| T-Cells | Tcell | 
| Excitatory Neurons | Excit |
| Inhibitory Neurons | Inhib |