Deconvolution Analysis
========

Louise A. Huuki-Myers
5/9/2022

### Details
[`bisque_deconvo.R`](bisque_deconvo.R): Run deconvolution analysis with Bisque 
(Jew et al., Nature Communications, 2013)  

* Used 5 region "pan" brain data from Tran, Maynard et al., Neuron, 2021 as sn-RNAseq reference  
* ten "broad" cell types  
* top 25 "mean ratio" marker genes (except 23 for Macro), recorded in `processed-data/02_deconvolution/PTSD_deconvolution_markers.csv`  
* results recorded in `processed-data/02_deconvolution/PTSD_deconvolution_results.csv` & `processed-data/02_deconvolution/est_prop_Bisque.Rdata`  

