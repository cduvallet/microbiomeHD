This directory contains version-controlled information about the data in [MicrobiomeHD](https://zenodo.org/record/569601#.WYMcMdPytTY).

# Files

## User-provided
* `zenodo_text.md` has the Zenodo description written in Markdown
* `dataset_info.yaml` has the dataset information. Note that this isn't necessarily the same as `../data/user_input/results_folder.yaml` because MicrobiomeHD has more datasets than was used in this analysis. Each dataset should have the following attributes:
    * `folder` - stem of the `*.tar.gz` file name on Zenodo
    * `paper` - DOI link to the original publication
    * `sample_size` - sample sizes

## Output files

`Makefile` makes `info_table.md` out of the information in `dataset_info.yaml` and the files uploaded to Zenodo.
Then it converts both `info_table.md` and `zenodo_text.md` into html files and concatenates them into `zenodo_description_final.html`. This is the description that you should copy and paste into the Zenodo website.

# Adding new datasets

To update Zenodo with a new dataset, upload the data to MicrobiomeHD and edit the `datasets_info.yaml` file in this folder accordingly.

Then, run `make` and copy the updated `zenodo_description_final.html` into Zenodo.
