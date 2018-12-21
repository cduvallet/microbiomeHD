# MicrobiomeHD: the human gut microbiome in Health and Disease

This repo contains the code to reproduce all of the analyses in "Meta analysis of microbiome studies identifies shared and disease-specific patterns", Duvallet et al. 2017.

The raw data is available on Zenodo: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.840333.svg)](https://doi.org/10.5281/zenodo.840333)

More information on the raw data on Zenodo is in the `db/` folder of this repo.

# Reproducing analyses

The paper is accompanied by a Makefile, which you can use to re-make all of the analyses, figures, and tables in the paper.

The supplementary files are included in the `final/supp-files/` folder of this
repo. Most of the folders in this repo are currently empty, and will get
populated with the results from `make`. The folders `data/user_input` and
`data/lit_search` are the other two folders with provided information.

`make` will download the data from Zenodo, clean and process the raw data,
perform all of the analyses in the paper, and make all of the figures and
tables.

Note that `make` does not do any of the random forest parameter search-related
analyses or figures. This part takes a very long time, and should be done as
a background process.

Finally, re-building the PhyloT tree takes some time and re-orders the genera
differently than in the paper (because there are multiple representations for
each tree). The tree used in the paper is provided in `data/analysis_results/`.
If you want to skip re-making the tree, you should run `make tree --touch`
before running `make`.

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
* `tree`: the phyloT tree used to order genera in final figures.

## Installing

To re-make all of the analyses, you'll first need to install the required
modules.

You should probably do this in a Python 2 virtual environment. Unfortunately, many of the packages I used are no longer available in conda, so if you use anaconda you'll need to first create an empty conda environment, [install pip](https://github.com/ContinuumIO/anaconda-issues/issues/1429), and then pip install the packages. If you don't use anaconda and/or have an alternative preferred way of making virtual environments, that's fine too (but I can't confirm that it will all work out). From the main directory, type:

```
conda create -n microbiomeHD python=2.7
source activate microbiomeHD
conda install pip
pip install -r requirements.txt
```

Note that all of these scripts were written in and for Python 2. Also, there have been many backward incompatible changes in some important modules used throughout, so you should install the old versions specified in the `requirements.txt` or else be plagued by many import errors.

You also need to install the NCBI EDirect command line tools for making
the tree. Instructions on how to do that are on the NCBI
[documentation](https://www.ncbi.nlm.nih.gov/books/NBK179288/).

Then you just run make:

`make`

And voila! A paper!

# Directory structure

This repo's structure follows what's recommended by [Cookie Cutter Data Science](https://drivendata.github.io/cookiecutter-data-science/).
Some of the files that are made and which I think might be useful to you
are already included in this repo. Other files are made by various scripts
in `src/`.

#### data

All data-related files are (or will be) in `data/`:

* `user_input`: user-inputted files, including:
  * `results_folder.yaml` file containing metadata on all the datasets [included]
  * `list_of_tar_files.txt` file that's used to download the raw data from Zenodo [included]
* `lit_search`: manual curation of the results reported in the original publications [included]
* `analysis_results`: files created by the analyses (e.g. q-values, random forest AUCs, etc) [made, except for the phyloT tree which is included]
* `raw_otu_tables`: the raw OTU tables as downloaded from Zenodo
* `clean_tables`: OTU tables and metadata in feather format, with "cleaned"
data (i.e. only samples with both metadata and 16S, OTUs and samples
with too few reads removed, etc)
* `tree`: files associated with the phyloT tree. Note that the final tree
(and its direct prerequisites) are included in this repo. If you want to
re-make the tree from scratch, delete any of these files before running
`make`. Making the tree is dicier and involves a manual step for you at
`http://phylot.biobyte.de/`. Re-making the tree will also change the order
of genera so that they no longer match the ordering in the paper exactly
(since the linear order of phylogenetic groups doesn't matter).

#### source code

All of the code is in the `src/` folder:

* `analysis`: all of the code used to perform any analyses
* `data`: data-related code, i.e. to download the raw data from Zenodo and to clean up the raw OTU tables and metadata files
* `final`: code used to make the final figures, tables, and supplementary files
* `util`: various functions and modules used in other scripts

#### figures, tables, and supplementary files

The Supplementary Files, Figures, and Tables are in the `final/` folder.
Some are made by `make`, others are included with this repo.

* `figures`: figures [made]
* `tables`: tables, in Markdown, tex, and tab-delimited formats [made]
* `supp-files`: supplementary files [included]
