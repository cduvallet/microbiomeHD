**Overview**

MicrobiomeHD is a standardized database of human gut microbiome studies in health and disease. This database includes publicly available 16S data from published case-control studies and their associated patient metadata. Raw sequencing data for each study was downloaded and processed through a standardized pipeline.

To be included in MicrobiomeHD, datasets have:

* publicly available raw sequencing data (fastq or fasta)
* publicly available metadata with at least case and control labels for each patient

Currently, MicrobiomeHD is focused on stool samples. Additional samples may be included in certain datasets, as indicated in the metadata.

**Files**

Additional information about the datasets included in this MicrobiomeHD release are in the MicrobiomeHD github repo [https://github.com/cduvallet/microbiomeHD](https://github.com/cduvallet/microbiomeHD), in the file *db/dataset_info.yaml*. Top-level identifiers correspond to dataset IDs labeled by disease_first-author.
For the most part, sample sizes in the yaml file are those that were described in the papers, and may not exactly reflect the actual data (due to missing/extra data, samples which didn't pass quality control, etc).

Each dataset was downloaded and processed through a standardized pipeline. The raw processing results are available in the \*.tar.gz files here. Each file has the same directory structure and files, as described in the pipeline documentation: [http://amplicon-sequencing-pipeline.readthedocs.io/en/latest/output.html](http://amplicon-sequencing-pipeline.readthedocs.io/en/latest/output.html).

Specific files of interest in each \*.tar.gz folder include:

* **summary_file.txt**: this file contains a summary of all parameters used to process the data
* **datasetID.metadata.txt**: the metadata associated with the samples. Note that some samples in the metadata may not have sequencing data, and vice versa.
* **RDP/datasetID.otu_table.100.denovo.rdp_assigned**: the 100% OTU tables with Latin taxonomic names assigned using the RDP classifier (c = 0.5).
* **datasetID.otu_seqs.100.fasta**: representative sequences for each OTU in the 100% OTU table. OTU labels in the OTU table end with `d__denovoID` - these denovoIDs correspond to the sequences in this file.
* **README.txt**: additional information about steps taken to download and process each dataset, as needed.

The raw data was acquired as described in the supplementary materials of Duvallet et al.'s "Meta analysis of microbiome studies identifies shared and disease-specific patterns" and, when available, the respective dataset README files.

Raw sequencing data was processed with the Alm lab's in-house 16S processing pipeline: [https://github.com/thomasgurry/amplicon_sequencing_pipeline](https://github.com/thomasgurry/amplicon_sequencing_pipeline)

Pipeline documentation is available at: [http://amplicon-sequencing-pipeline.readthedocs.io/](http://amplicon-sequencing-pipeline.readthedocs.io/)

Metadata was extracted from the original papers and/or data sources, and formatted manually. When possible, these steps are documented in each dataset's
associated README.txt file.

**Contributing**

MicrobiomeHD is a resource that can be used to extract disease-specific microbiome signals in individual case-control studies.
Many microbes respond non-specifically to health and disease, and the majority of bacterial associations within individual studies overlap with this non-specific response.
Researchers should cross-check their results with the data presented here to ensure that their identified microbial associations are specific to their disease under study.

We provide an updated list of non-specific microbes here, as well as the raw OTU tables for anyone who wishes to reproduce and adapt this analysis to their study question.

If you would like to include your case-control dataset in MicrobiomeHD, please email ejalm[at]mit.edu and duvallet[at]mit.edu.

For us to process your data through our standard pipeline, you will need to provide the following files and information about your data:

* raw sequencing data in fastq or fasta format (preferably fastq)
* information about which processing steps will be required (e.g. removing primers or barcodes, merging paired-end reads, etc)
* sample IDs associated with the sequencing data (either mapped to barcodes still in the sequences, or to each de-multiplexed sequencing file)
* case/control metadata of each sample
* other relevant metadata (e.g. sampling site, if not all samples are stool; sampling time point, if multiple samples per patient were taken; etc)

By using MicrobiomeHD in your own analyses, you agree to contribute your dataset to this database and to make your raw sequencing data (i.e. fastq files) publicly available.

**Citing MicrobiomeHD**

The MicrobiomeHD database and original publications for each of these datasets are described in Duvallet et al. (2017): [http://dx.doi.org/10.1038/s41467-017-01973-8](http://dx.doi.org/10.1038/s41467-017-01973-8)

Duvallet, C., Gibbons, S. M., Gurry, T., Irizarry, R. A., & Alm, E. J. (2017).
Meta-analysis of gut microbiome studies identifies disease-specific and shared
responses. *Nature communications*, 8(1), 1784.

If you use any of these datasets in your analysis, please cite both MicrobiomeHD (Duvallet et al. (2017)) and the original publication for each dataset that you use.

The code used to process and analyze this data in the paper is available on github: [https://github.com/cduvallet/microbiomeHD](https://github.com/cduvallet/microbiomeHD)

**Files**

*Data files*

**file-S3.nonspecific_genera.txt**: Supplemental Table 3 from Duvallet et al. (2017), listing the non-specific health- and disease-associated microbes.   
**dataset_info.yaml**: yaml file with additional dataset metadata.   

*Datasets*

Note that MicrobiomeHD contains all 28 datasets from Duvallet et al. (2017), as well as additional datasets which did not meet the inclusion criteria for the meta-analysis presented in the paper. Additional information about the datasets included in this MicrobiomeHD release are in the original publications and the MicrobiomeHD github repo https://github.com/cduvallet/microbiomeHD, and in the file *dataset_info.yaml*.

The sample sizes listed here reflect what was reported in the original publications. Some may have discrepancies between what is reported and what is in the actual data due to missing data, quality issues, barcode mismatches, etc.
