### Make file to reproduce data and figures from "Universal patterns of gut
### microbiome shifts in disease" by Duvallet et al.

## User notes
# Must have list_of_tar_files.txt and results_folders.yaml files
# datasetIDs cannot have slashes or periods in them
# If you change the raw data, clean OTU tables, or clean metadata files,
# make sure to delete both the clean OTU and metadata tables

## Some Makefile notes
# Automatic variables: https://www.gnu.org/software/make/manual/html_node/Automatic-Variables.html#Automatic-Variables
# $@ is the target file
# $* is the stem that files have in common




all: data analysis #all_figures

### User inputs
tar_files = data/list_of_tar_files.txt
yaml_file = data/results_folders.yaml

### Define target data files
# define the tar files from the list_of_tar_files.txt file
raw_tar_files := $(shell cat data/list_of_tar_files.txt)
# define the clean OTU table names from the dataset IDs in the results_folders.yaml
clean_otu_tables := $(shell grep -v '^    ' $(yaml_file) | grep -v '^\#' | sed 's/:/.otu_table.clean/g' | sed 's/^/data\/clean_tables\//g')
# define the metadata file names from the dataset IDs in the results_folders.yaml
clean_metadata_files := $(shell grep -v '^    ' $(yaml_file) | grep -v '^\#' | sed 's/:/.metadata.clean/g' | sed 's/^/data\/clean_tables\//g')

dataset_info = data/datasets_info/datasets_info.txt
proc_info = data/datasets_info/datasets_info.processing.txt
manual_meta_analysis = data/lit_search/literature_based_meta_analysis.txt

raw_data: $(raw_tar_files)
clean_data: $(clean_otu_tables) #$(clean_metadata_files)
data_info: $(dataset_info) $(manual_meta_analysis)

data: raw_data clean_data data_info

## 1. Pull raw OTU tables from somewhere (prob Zenodo? for now, my other folder).
# output: tar files in a directory called data/raw_otu_tables,
$(raw_tar_files): src/data/copy_tar_folders.sh
	src/data/copy_tar_folders.sh $@

## 2. Clean the raw OTU tables and metadata files
# Note: technically these clean files should depend on the raw_tar_files,
# but because the stem of the files doesn't necessarily match, including
# raw_tar_files as a prerequisite means that make can no longer parallelize
# the cleaning steps. So if you changed the raw data for one of the OTU tables,
# you need to delete it so that make knows to re-process it.
# Note: If my results_folders were better labeled, I could simply
# write a rule like %.otu_table.clean : clean.py raw_data/%.tar.gz

# This code cleans both the OTU and metadata files,
# and writes both *.otu_table.clean and *.metadata.clean
$(clean_otu_tables): src/data/clean_otu_and_metadata.py $(yaml_file)
	python src/data/clean_otu_and_metadata.py data/raw_otu_tables $(yaml_file) $@

## 3. Get info about datasets
# input: files in the clean_data/ folder, and a script to get info for all datasets
# output: datasets_info.txt
$(dataset_info): src/data/dataset_info.py
	python src/data/dataset_info.py $(yaml_file) data/raw_otu_tables data/clean_tables $(dataset_info) $(proc_info)

$(proc_info): src/data/dataset_info.py
		python src/data/dataset_info.py $(yaml_file) data/raw_otu_tables data/clean_tables $(dataset_info) $(proc_info)

## 4. Manual meta-analysis
# Maybe pull this from the internet? Maybe just have it echo the link to download it from? Or, like, "grab it from the paper"?
$(manual_meta_analysis):
	echo -e "Download this from the supplement or something"

## Define the target analysis files

qvalues = data/analysis_results/q-val_all_results.mean.kruskal-wallis.case-control.txt
meta_qvalues = data/analysis_results/meta.counting.q-0.05.disease_wise.txt \
               data/analysis_results/meta.counting.q-0.05.2_diseases.across_all_diseases.txt

analysis: $(qvalues) $(meta_qvalues)

## 1. q-values files for all genera across all studies
$(qvalues): src/analysis/get_qvalues.py $(clean_otu_tables) $(clean_metadata_files)
	python src/analysis/get_qvalues.py data/clean_tables data/analysis_results

## 2. meta-analysis results (+/- 1's) for each disease
$(meta_qvalues): src/analysis/meta_analyze.py $(qvalues)
	python src/analysis/meta_analyze.py $(qvalues) data/analysis_results 0.05 2

## 3. overall meta-analysis results (maybe combined with 2? I don't remember how I structured the codes and the files)

## 4. phyloT stuff
# this outputs a file, right?

## 5. alpha diversities

## 6. random forest results

## 7. random forest parameter search

### make figures
# Just go through the figures in the directory structure business

# Also make the tables from datasets_info.py

### make paper
# Just grab from my latest repo? Dunno.
