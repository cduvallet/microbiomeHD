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




all: data analysis tree #all_figures

### User inputs
tar_files = data/user_input/list_of_tar_files.txt
yaml_file = data/user_input/results_folders.yaml

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
clean_data: $(clean_otu_tables) $(clean_metadata_files)
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

$(clean_metadata_files): $(clean_otu_tables)
	## Recover from the removal of $@
	## i.e. if the metadata file is deleted but the OTU table still is unchanged
	@if test -f $@; then :; else \
		rm -f $<; \
		$(MAKE) $(AM_MAKEFLAGS) $<; \
  	fi

## 3. Get info about datasets
# input: files in the clean_data/ folder, and a script to get info for all datasets
# output: datasets_info.txt
$(dataset_info): src/data/dataset_info.py $(yaml_file) $(clean_otu_tables) $(clean_metadata_files)
	python src/data/dataset_info.py $(yaml_file) data/raw_otu_tables data/clean_tables $(dataset_info) $(proc_info)

$(proc_info): $(dataset_info)
	## Recover from the removal of $@
	## i.e. if the proc_info file is deleted but the dataset_info is unchanged,
	## since dataset_info.py makes both
	@if test -f $@; then :; else \
		rm -f $(dataset_info); \
		$(MAKE) $(AM_MAKEFLAGS) $(dataset_info); \
  	fi

## 4. Manual meta-analysis
# Maybe pull this from the internet? Maybe just have it echo the link to download it from? Or, like, "grab it from the paper"?
$(manual_meta_analysis):
	echo -e "Download this from the supplement or something"

## Define the target analysis files
qvalues = data/analysis_results/q-val_all_results.mean.kruskal-wallis.case-control.txt
meta_qvalues = data/analysis_results/meta.counting.q-0.05.disease_wise.txt
overall_qvalues = data/analysis_results/meta.counting.q-0.05.2_diseases.across_all_diseases.txt
phyloT_tree = data/analysis_results/genus_tree.tre


analysis: $(qvalues) $(meta_qvalues) $(overall_qvalues)

## 1. q-values files for all genera across all studies
$(qvalues): src/analysis/get_qvalues.py $(clean_otu_tables) $(clean_metadata_files)
	python src/analysis/get_qvalues.py data/clean_tables data/analysis_results

## 2. meta-analysis results (+/- 1's) for each disease
$(meta_qvalues): src/analysis/meta_analyze.py $(qvalues)
	python src/analysis/meta_analyze.py $(qvalues) data/analysis_results 0.05 2

## 3. overall meta-analysis results, made at the same time as meta_qvalues
$(overall_qvalues): $(meta_qvalues)
	@if test -f $@; then :; else \
		rm -f $(meta_qvalues); \
		$(MAKE) $(AM_MAKEFLAGS) $(meta_qvalues); \
	fi

## 5. alpha diversities

## 6. random forest results

## 7. random forest parameter search

## Tree stuff
genera_file = data/analysis_results/genera.tmp
ncbi_file = data/analysis_results/ncbi_ids.tmp
clean_ncbi = data/analysis_results/ncbi_ids.clean.for_phyloT
phyloT_file = data/user_input/phyloT_tree.newick
final_tree_file = data/analysis_results/phyloT_tree.updated.newick

tree: $(final_tree_file)

# Grab genus names from the qvalues file
$(genera_file): src/analysis/genera_from_qvalues.py $(qvalues)
	src/analysis/genera_from_qvalues.py $(qvalues) $(genera_file)

# Get NCBI IDs using esearch
$(ncbi_file): src/analysis/get_ncbi_IDs.sh $(genera_file)
	./src/analysis/get_ncbi_IDs.sh $(genera_file) $(ncbi_file)

# Clean up the NCBI IDs (i.e. remove non-Bacteria things)
# and keep just the genus IDs
$(clean_ncbi): src/analysis/clean_ncbi.py $(ncbi_file)
	src/analysis/clean_ncbi.py $(ncbi_file) data/analysis_results/ncbi_ids.clean.tmp $(clean_ncbi)

# Manual step: go to the phyloT website and make the tree
$(phyloT_file): $(clean_ncbi)
	read -n1 -p "Go to http://phylot.biobyte.de/ and generate tree from ${clean_ncbi}. Press any key to continue once you've added the tree file in ${phyloT_file}. "

# Manually edit tree with genera which didn't have NCBI IDs
$(final_tree_file): src/analysis/update_tree.py $(phyloT_file) $(genera_file)
	src/analysis/update_tree.py $(genera_file) $(phyloT_file) $(final_tree_file)

# Re-order the qvalues and meta-analysis files phylogenetically


### make figures
# Just go through the figures in the directory structure business

# Also make the tables from datasets_info.py

### make paper
# Just grab from my latest repo? Dunno.
