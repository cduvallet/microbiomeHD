# MicrobiomeHD: the human gut microbiome in Health and Disease

This repo contains the code to reproduce all of the analyses in "Meta analysis of microbiome studies identifies shared and disease-specific patterns", Duvallet et al. 2017.

The raw data is available on Zenodo: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.569601.svg)](https://doi.org/10.5281/zenodo.569601)


# Reproducing analyses

[in progress]

The paper is accompanied by a Makefile, which you can use to re-make all of the analyses, figures, and tables in the paper.

`make` will download the data from Zenodo, clean and process the raw data, perform all of the analyses in the paper, and
make all of the figures and tables.

Other things you can `make`:

* `figures`: all of the figures
  * `main_figures`: just the figures in the main text
  * `supp_figures`: supplementary figures
* `tables`: all of the tables, in both Markdown and tex formats
* `analysis`: all of the analysis files, but none of the figures or tables
* `rf_params`: the random forest parameter search analysis. Note that this takes a very long time, and should be done
as a background process. It is not included in any of the other `make` commands.


# Directory structure

This repo's structure follows what's recommended by [Cookie Cutter Data Science](https://drivendata.github.io/cookiecutter-data-science/).
Some of the files have been included in this repo, others are made by various scripts.

All data-related files are in `data/`:

* user_input: user-inputted files, including:
  * `results_folder.yaml` file containing metadata on all the datasets [included]
  * `list_of_tar_files.txt` file that's used to download the raw data from Zenodo [included]
* lit_search: manual curation of the results reported in the original publications [included]
* analysis_results: files created by the analyses (e.g. q-values, random forest AUCs, etc) [made]

All of the code is in the `scr/` folder:

* `analysis`: all of the code used to perform any analyses
* `data`: data-related code, i.e. to download the raw data from Zenodo and to clean up the raw OTU tables and metadata files
* `final`: code used to make the final figures, tables, and supplementary files
* `util`: various functions and modules used in other scripts

The Supplementary Files, Figures, and Tables are in the `final/` folder. The Supplementary Files are provided in this repo;
the Tables and Figures are made by `make`.
