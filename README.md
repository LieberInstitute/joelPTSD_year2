[![DOI](https://zenodo.org/badge/373568668.svg)](https://zenodo.org/badge/latestdoi/373568668)

# joelPTSD_year2

This reposotory contains the custom code for processing the bulk RNA sequencing (RNA-seq) data from US National Institutes of Health grant [R01MH117291](https://reporter.nih.gov/search/B68mnmkOgkGdPZcs09B2zw/project-details/10400227) whose LIBD principal investigator is [Joel E. Kleinman](https://www.libd.org/team/joel-kleinman).

If you have any questions about the code on this GitHub repository, please ask questions publicly at https://github.com/LieberInstitute/joelPTSD_year2/issues such that everyone can benefit from the answers. Thank you!

For citing this work, use the following information in [BibTeX](http://www.bibtex.org/) format:

```
## TODO: add citation info
```

# Main files

## RNA-seq processing


The [RNAseq-pipeline](https://github.com/LieberInstitute/RNAseq-pipeline) (DOI: [10.5281/zenodo.7348660](https://doi.org/10.5281/zenodo.7348660)) was used to preprocess the first batch of sequencing reads. The later batches of data wasere processed with [SPEAQeasy](https://github.com/LieberInstitute/SPEAQeasy) (PMID [33932985](https://pubmed.ncbi.nlm.nih.gov/33932985/)), the continuation of the RNAseq-pipeline project. Both workflows align the sequencing reads to the human genome sequence GRCh38  using the [HISAT2 aligner](http://daehwankimlab.github.io/hisat2/) (PMID [25751142](https://pubmed.ncbi.nlm.nih.gov/25751142/)). The pipeline produced [SummarizedExperiment](https://bioconductor.org/packages/SummarizedExperiment/) R objects (PMID [25633503](https://pubmed.ncbi.nlm.nih.gov/25633503/)) with quantification of all genomic features as annotated by [GENCODE release 25](https://www.gencodegenes.org/human/release_25.html) (PMID [33270111](https://pubmed.ncbi.nlm.nih.gov/33270111/)). The pipeline also provided coverage data for [unannotated expressed regions](https://bioconductor.org/packages/derfinder/) (PMID [27694310](https://pubmed.ncbi.nlm.nih.gov/27694310/)) and all exon-exon junctions found in the RNA-seq read alignments. Custom R code (from this repository) was used to merge data across the batches processed with the RNAseq-pipeline and SPEAQeasy, forming a single set of SummarizedExperiment R objects for downstream analysis.

### Main scripts

TODO by @nick-eagles

## Cell type deconvolution

TODO by @lahuuki


### Main scripts

TODO by @lahuuki

## Internal

JHPCE location: `/dcs04/lieber/lcolladotor/joelPTSD_LIBD3040/joelPTSD_year2`

