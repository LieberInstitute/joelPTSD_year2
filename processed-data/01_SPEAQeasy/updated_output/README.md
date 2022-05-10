Data contained here, by directory:

- `processed_data`: final processed data for year 2. Because of the difference in processing date with year 1, year 2 was processed with a more recent version of the computational pipeline now called [SPEAQeasy](https://github.com/LieberInstitute/SPEAQeasy). The `Rdata` files contain `RangedSummarizedExperiment` objects for each feature type (genes, exons, transcripts, and exon-exon junctions) described in more detail [here](http://research.libd.org/SPEAQeasy/outputs.html#main-outputs). A CSV table of collected metrics and other stats for each sample is also here, which contains the same information as the `colData` of each R object.

- `processed_data/merged_with_Y1`: final processed R objects containing all samples across year 1 and 2.

- `fastq`: raw sequencer outputs (FASTQ format) for year 2.


From Geo Pertea:

* Year2 dataset with demo data added: `/dcs04/lieber/lcolladotor/joelPTSD_LIBD3040/joelPTSD_year2/processed-data/01_SPEAQeasy/updated_output/year_2/rse_*.y2_n238.Rdata`
* Year1+Year2 merged RSEs: `/dcs04/lieber/lcolladotor/joelPTSD_LIBD3040/joelPTSD_year2/processed-data/01_SPEAQeasy/updated_output/merged/rse_*.y1y2_n463.Rdata`.
* Note that the Y1 data in this merged file was already included in the large merged RSEs as part of the "PTSD_Brainomics" dataset.
* The pheno data from Ran that I used for this Year 2 batch is here: `/dcs04/lieber/lcolladotor/joelPTSD_LIBD3040/joelPTSD_year2/pheno_data/Joel_PTSD_Y3_RNA_Samples_Info.csv`
