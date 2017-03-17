### Make file to reproduce data and figures from "Universal patterns of gut
### microbiome shifts in disease" by Duvallet et al.

## User notes
# Must have list_of_tar_files.txt and results_folders.yaml files
# datasetIDs cannot have slashes or periods in them

## Some Makefile notes
# Automatic variables: https://www.gnu.org/software/make/manual/html_node/Automatic-Variables.html#Automatic-Variables
# $@ is the target file
# $* is the stem that files have in common


all: data #all_figures

## make data

## File definitions
tar_files = data/list_of_tar_files.txt
yaml_file = data/results_folders.yaml

# define the tar files from the list_of_tar_files.txt file
raw_tar_files := $(shell cat data/list_of_tar_files.txt)
# define the clean OTU table names from the dataset IDs in the results_folders.yaml
clean_otu_tables := $(shell grep -v '^    ' $(yaml_file) | grep -v '^\#' | sed 's/:/.otu_table.clean/g' | sed 's/^/data\/clean_tables\//g')
# define the metadata file names from the dataset IDs in the results_folders.yaml
clean_metadata_files := $(shell grep -v '^    ' $(yaml_file) | grep -v '^\#' | sed 's/:/.metadata.clean/g' | sed 's/^/data\/clean_tables\//g')
# clean_otu_tables = something # $(call a python function to make list of clean tables that we expect)
# clean_metadata = something_similar_to_above # data/clean_tables/*.metadata.clean
# dataset_info = data/clean_tables/datasets_info.txt
# analysis_results = manually_defined_files_here

data: $(raw_tar_files) $(clean_otu_tables) $(clean_metadata_files)

## 1. Pull raw OTU tables from somewhere (prob Zenodo? for now, my other folder).
# output: tar files in a directory called data/raw_otu_tables,
$(raw_tar_files): src/data/copy_tar_folders.sh
	src/data/copy_tar_folders.sh $@

## 2. Clean the raw OTU tables and metadata files
# Targets: clean otu table file names
# Rule: python clean_otu_tables.py raw_data_dir results_folder.yaml clean_otu_fname
$(clean_otu_tables): src/data/clean_otu_tables.py
	python src/data/clean_otu_tables.py data/raw_otu_tables data/results_folders.yaml $@

# Targets: clean metadata file names
# Rule: python clean_metadata.py raw_data_dir results_folder.yaml clean_metadata_fname
# Note: this also updates the OTU tables!
$(clean_metadata_files): src/data/clean_metadata.py $(clean_otu_tables)
	python src/data/clean_metadata.py data/raw_otu_tables data/results_folders.yaml $@
# input: folders in the data/raw_otu_tables directory, and a script to clean one dataset at a time (given a yaml file with the location of the metadata, etc?)
# output: *.otu_table.clean and *.metadata.clean files with the same prefix as the folder prefixes, in a directory called data/clean_data


## 3. Get info about datasets
# input: files in the clean_data/ folder, and a script to get info for all datasets
# output: datasets_info.txt, maybe some latex strings too? (nah)

## 4. Manual meta-analysis
# Maybe pull this from the internet? Maybe just have it echo the link to download it from? Or, like, "grab it from the paper"?

## make analyses

## 1. q-values files for all genera across all studies

## 2. meta-analysis results (+/- 1's) for each disease

## 3. overall meta-analysis results (maybe combined with 2? I don't remember how I structured the codes and the files)

## 4. phyloT stuff
# this outputs a file, right?

## 5. alpha diversities


### make figures
# Just go through the figures in the directory structure business

### make paper
# Just grab from my latest repo? Dunno.
