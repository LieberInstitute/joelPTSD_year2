Data contained here, by directory:

- `processed_data`: final processed data for year 2. Because of the difference in processing date with year 1, year 2 was processed with a more recent version of the computational pipeline now called [SPEAQeasy](https://github.com/LieberInstitute/SPEAQeasy). The `Rdata` files contain `RangedSummarizedExperiment` objects for each feature type (genes, exons, transcripts, and exon-exon junctions) described in more detail [here](http://research.libd.org/SPEAQeasy/outputs.html#main-outputs). A CSV table of collected metrics and other stats for each sample is also here, which contains the same information as the `colData` of each R object.

- `processed_data/merged_with_Y1`: final processed R objects containing all samples across year 1 and 2.

- `fastq`: raw sequencer outputs (FASTQ format) for year 2.
