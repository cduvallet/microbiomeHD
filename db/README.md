This directory contains version-controlled information about the data in [MicrobiomeHD](https://doi.org/10.5281/zenodo.797943).

`make` just compiles a few different sources of information to update the
description that's on Zenodo.

# Files

## User-provided

*Notes for contributors and MicrobiomeHD maintainer:*

* `dataset_info.yaml` has the dataset information. Note that this isn't necessarily the same as `../data/user_input/results_folder.yaml` because MicrobiomeHD has more datasets than was used in this analysis. Each dataset should have the following attributes:
    * `folder` - stem of the `*.tar.gz` file name on Zenodo
    * `paper` - DOI link to the original publication
    * `sample_size` - sample sizes
* `zenodo_text.md` has the Zenodo description written in Markdown

## Output files

*Notes for MicrobiomeHD maintainer:*

`Makefile` makes `info_table.md` out of the information in `dataset_info.yaml` and the files uploaded to Zenodo.
Then it converts both `info_table.md` and `zenodo_text.md` into html files and concatenates them into `zenodo_description_final.html`. This is the description that is copied for my Zenodo description.

# Adding new datasets

*Notes for contributors:*

If you want to include your dataset in MicrobiomeHD, email Claire at
duvallet[at]mit[dot]edu and Eric at ejalm[at]mit[dot]edu. We'll work with
you to process and upload your data so that it can be compared with the other
datasets in MicrobiomeHD.

*Notes for MicrobiomeHD maintainer:*

To then update Zenodo with a new dataset, upload the data to MicrobiomeHD and edit the `datasets_info.yaml` file in this folder accordingly.

Then, run `make` and copy the updated `zenodo_description_final.html` into Zenodo.
