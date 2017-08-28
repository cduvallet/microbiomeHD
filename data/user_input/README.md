# User-inputted data files

#### results_folders.yaml

This file contains dataset-related metadata, and provides the interface between the raw
OTU tables and all downstream processing steps.

Top-level keys are the datasetID as it will be used in all downstream code and analyses.

Each dataset should have the following keys:

* `folder`: <folder>.tar.gz is the name of the file that is downloaded from Zenodo.
When extracted, this folder should contain `dataset.metadata.txt` file and
`RDP/dataset.otu_table.100.denovo.rdp_assigned` files. The value for `folder` is of the
format: `<dataset>_results`, where `dataset` is the stem file name for all files in the
results folder. This value of `dataset` will be converted to the top-level datasetID when
data is cleaned.
* `data_source`, `metadata_source`: information about where the raw sequencing data and
metadata were acquired from. This information is extracted to make Supplementary Table 4.
* `region`, `sequencer`, `year`: dataset metadata used to populate Table 1 and Supplementary
Table 2.

Other keys are optional and provide helpful notes about each dataset.

#### list_of_tar_files.txt

The list of tar files to download from Zenodo. This file should have no header, and have each
file name on a new line.

```
>claire:~/github/microbiomeHD/data/user_input$ cat list_of_tar_files.txt | head -n 4
cdi_youngster_results.tar.gz
hiv_noguerajulian_results.tar.gz
hiv_lozupone_results.tar.gz
t1d_alkanani_results.tar.gz
```

#### split_cases_datasets.txt

This file contains the dataset IDs for which separate analyses should be done
for each case group (i.e. in the Supplementary Discussion "Stratifying
 heterogenous case groups shows consistent disease-specific signals").

 It should have each dataset ID on its own line.
