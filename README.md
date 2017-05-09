# MicrobiomeHD: the human gut microbiome in Health and Disease

This repo contains the code to reproduce all of the analyses in "Meta analysis of microbiome studies identifies shared and disease-specific patterns", Duvallet et al. 2017.

The raw data is available on Zenodo: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.569601.svg)](https://doi.org/10.5281/zenodo.569601)

More information on the raw data on Zenodo is in `data/user_input/results_folders.yaml`.

# Reproducing analyses

The paper is accompanied by a Makefile, which you can use to re-make all of the analyses, figures, and tables in the paper.

`make` will download the data from Zenodo, clean and process the raw data, perform all of the analyses in the paper, and
make all of the figures and tables.

Note that `make` does not do any of the random forest parameter search-related
analyses or figures. This part takes a very long time, and should be done as
a background process.

Other things you can `make` separately:

* `figures`: all of the figures
  * `main_figures`: just the figures in the main text
  * `supp_figures`: supplementary figures
* `tables`: all of the tables, in both Markdown and tex formats
* `analysis`: all of the analysis files, but none of the figures or tables
* `rf_params`: the random forest parameter search analysis. Note that this
is not included in any of the other `make` commands.
* `supp_files`: the supplementary files, which are also included in the repo
and don't technically need to be re-made

To re-make all of the analyses, you'll first need to install the required
modules. From this main directory, you can type:

`pip install -r requirements.txt`

Then you just run make:

`make`

And voila! A paper!

# Directory structure

This repo's structure follows what's recommended by [Cookie Cutter Data Science](https://drivendata.github.io/cookiecutter-data-science/).
Some of the files that are made and which I think might be useful to you
are already included in this repo. Other files are made by various scripts
in `src/`.

#### data

All data-related files are in `data/`:

* `user_input`: user-inputted files, including:
  * `results_folder.yaml` file containing metadata on all the datasets [included]
  * `list_of_tar_files.txt` file that's used to download the raw data from Zenodo [included]
* `lit_search`: manual curation of the results reported in the original publications [included]
* `analysis_results`: files created by the analyses (e.g. q-values, random forest AUCs, etc) [made, except for the phyloT tree which is included]

#### source code

All of the code is in the `src/` folder:

* `analysis`: all of the code used to perform any analyses
* `data`: data-related code, i.e. to download the raw data from Zenodo and to clean up the raw OTU tables and metadata files
* `final`: code used to make the final figures, tables, and supplementary files
* `util`: various functions and modules used in other scripts

#### figures, tables, and supplementary files

The Supplementary Files, Figures, and Tables are in the `final/` folder.

* `figures`: figures [made]
* `tables`: tables, in Markdown, tex, and tab-delimited formats [made]
* `supp-files`: supplementary files [included]
