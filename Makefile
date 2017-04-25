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
# Multiple targets for one rule: https://www.gnu.org/software/automake/manual/html_node/Multiple-Outputs.html

## Some common files
util = src/util/util.py
fileio = src/util/FileIO.py
summaryparser = src/util/SummaryParser.py


all: data analysis tree #all_figures

##### DOWNLOAD AND CLEAN DATA #####

### User inputs
tar_files = data/user_input/list_of_tar_files.txt
yaml_file = data/user_input/results_folders.yaml

### Define target data files
# define the tar files from the list_of_tar_files.txt file
raw_tar_files := $(shell cat data/user_input/list_of_tar_files.txt)
# define the clean OTU table names from the dataset IDs in the results_folders.yaml
clean_otu_tables := $(shell grep -v '^    ' $(yaml_file) | grep -v '^\#' | sed 's/:/.otu_table.clean.feather/g' | sed 's/^/data\/clean_tables\//g')
# define the metadata file names from the dataset IDs in the results_folders.yaml
clean_metadata_files := $(shell grep -v '^    ' $(yaml_file) | grep -v '^\#' | sed 's/:/.metadata.clean.feather/g' | sed 's/^/data\/clean_tables\//g')

manual_meta_analysis = data/lit_search/literature_based_meta_analysis.txt

raw_data: $(raw_tar_files)
clean_data: $(clean_otu_tables) $(clean_metadata_files)

data: raw_data clean_data $(manual_meta_analysis)

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

## 3. Manual meta-analysis
# Maybe pull this from the internet? Maybe just have it echo the link to download it from? Or, like, "grab it from the paper"?
$(manual_meta_analysis):
	echo -e "Download this from the supplement or something"

##### ANALYSES #####

## Define the target analysis files
qvalues = data/analysis_results/q-val_all_results.mean.kruskal-wallis.case-control.txt
meta_qvalues = data/analysis_results/meta.counting.q-0.05.disease_wise.txt
overall_qvalues = data/analysis_results/meta.counting.q-0.05.2_diseases.across_all_diseases.txt

alpha_divs = data/analysis_results/alpha_diversity.txt
alpha_pvals = data/analysis_results/alpha_diversity.pvalues.txt

rf_results = data/analysis_results/rf_results.txt
rf_param_search = data/analysis_results/rf_results.parameter_search.txt

analysis: $(qvalues) $(meta_qvalues) $(overall_qvalues) $(alpha_divs) $(alpha_pvals) $(rf_results) $(rf_param_search)

## 1. q-values files for all genera across all studies
qvals: $(qvalues) $(meta_qvalues) $(overall_qvalues)
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

## 4. alpha diversities
alpha: $(alpha_divs) $(alpha_pvals)
$(alpha_divs): src/analysis/alpha_diversity.py $(clean_otu_tables) $(clean_metadata_files)
	python src/analysis/alpha_diversity.py data/clean_tables \
	$(alpha_divs) $(alpha_pvals)

$(alpha_pvals): $(alpha_divs)
	@if test -f $@; then :; else \
		rm -f $(alpha_pvals); \
		$(MAKE) $(AM_MAKEFLAGS) $(alpha_pvals); \
	fi

## 5. random forest results
rf: $(rf_results) $(rf_param_search)

$(rf_results): src/analysis/classifiers.py $(clean_otu_tables) $(clean_metadata_files) src/util/util.py
	python src/analysis/classifiers.py data/clean_tables \
	$(rf_results)

## 6. random forest parameter search
$(rf_param_search): src/analysis/classifiers_parameters.py $(clean_otu_tables) $(clean_metadata_files)
	python src/analysis/classifiers_parameters.py data/clean_tables $(rf_param_search)

##### PHYLOT TREE STUFF #####
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

##### PREP FOR PLOTTING #####

## Prepare q-value files for plotting
# Re-order rows in qvalues and meta-analysis files phylogenetically
# Keep only genera which were significant in at least one study
# Also make the logfold change qvalues file
qvalues_clean = $(subst txt,sig_ordered.txt,$(qvalues))
meta_clean = $(subst txt,sig_ordered.txt,$(meta_qvalues))
overall_clean = $(subst txt,sig_ordered.txt,$(overall_qvalues))
logfold = $(subst txt,log2change.sig_ordered.txt,$(qvalues))

for_plotting: $(qvalues_clean) $(meta_clean) $(overall_clean) $(logfold)

$(qvalues_clean): src/analysis/reorder_qvalues.py $(qvalues) $(meta_qvalues) $(overall_qvalues) $(final_tree_file)
	python src/analysis/reorder_qvalues.py $(qvalues) --qthresh 0.05  $(meta_qvalues) $(overall_qvalues) $(final_tree_file)

$(meta_clean): $(qvalues_clean)
	@if test -f $@; then :; else \
		rm -f $(meta_clean); \
		$(MAKE) $(AM_MAKEFLAGS) $(qvalues_clean); \
	fi

$(overall_clean): $(qvalues_clean)
	@if test -f $@; then :; else \
		rm -f $(overall_clean); \
		$(MAKE) $(AM_MAKEFLAGS) $(qvalues_clean); \
	fi

## Calculate logfold change for all of the "clean" genera
# (i.e sig in at least one study, phylogenetically ordered)
$(logfold): src/analysis/logfold_effect.py $(qvalues_clean) $(clean_otu_tables) $(clean_metadata_files)
	python src/analysis/logfold_effect.py data/clean_tables $(qvalues_clean) $(logfold)

##### FIGURES AND TABLES #####

## Tables
dataset_info = data/datasets_info/datasets_info.txt
# Table 1 (main text) has the datasets, controls, cases, and references
table1 = final/tables/table1.tex
# Table 2 (supplement) has the same things, plus some info about the sequencers
table2 = final/tables/table2.tex
# Table 3 has the processing parameters
table3 = final/tables/table3.tex
# Table 4 has the data and metadata sources
table4 = final/tables/table4.tex

tables: $(table1) $(table2) $(table3) $(table4)

$(table1): src/figures-tables/tables-1-2.datasets_info.py $(yaml_file) $(clean_otu_tables) $(clean_metadata_files)
	python src/figures-tables/tables-1-2.datasets_info.py $(yaml_file) data/raw_otu_tables data/clean_tables $(dataset_info) $(table1) $(table2)

$(table2): $(table1)
	@if test -f $@; then :; else \
		rm -f $(table1); \
		$(MAKE) $(AM_MAKEFLAGS) $(table1); \
	fi

$(table3): src/figures-tables/tables-3-4.processing_info.py $(yaml_file)
	python src/figures-tables/tables-3-4.processing_info.py $(yaml_file) data/raw_otu_tables $(table3) $(table4)

$(table4): $(table3)
	@if test -f $@; then :; else \
		rm -f $(table3); \
		$(MAKE) $(AM_MAKEFLAGS) $(table3); \
	fi

### make figures
# Just go through the figures in the directory structure business

# Also make the tables from datasets_info.py

### make paper
# Just grab from my latest repo? Dunno.
