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

all: data analysis figures tables supp_files

# All of the files associated with the RF parameter search
rf_params: rf_param_search rf_param_figures

### User inputs - these files are not made, they should already exist!
tar_files = data/user_input/list_of_tar_files.txt
yaml_file = data/user_input/results_folders.yaml
manual_meta_analysis = data/lit_search/literature_based_meta_analysis.txt

###################################
##### DOWNLOAD AND CLEAN DATA #####
###################################

### Define target data files
# define the tar files from the list_of_tar_files.txt file
raw_tar_files := $(addprefix data/raw_otu_tables/,$(shell cat data/user_input/list_of_tar_files.txt))
# define the clean OTU table names from the dataset IDs in the results_folders.yaml
clean_otu_tables := $(shell grep -v '^    ' $(yaml_file) | grep -v '^\#' | sed 's/:/.otu_table.clean.feather/g' | sed 's/^/data\/clean_tables\//g')
# define the metadata file names from the dataset IDs in the results_folders.yaml
clean_metadata_files := $(shell grep -v '^    ' $(yaml_file) | grep -v '^\#' | sed 's/:/.metadata.clean.feather/g' | sed 's/^/data\/clean_tables\//g')

raw_data: $(raw_tar_files)
clean_data: $(clean_otu_tables) $(clean_metadata_files)
data: raw_data clean_data $(manual_meta_analysis)

## 1. Download the raw tar.gz files from Zenodo into data/raw_otu_tables,
## only if the file doesn't already exist. Also extract the files.
# Note: when I download from Zenodo, the file date corresponds to the day
# I uploaded the data to Zenodo (May 3). Need to touch the file to update the
# modified date so that make doesn't re-make these files all the time.
$(raw_tar_files): data/user_input/list_of_tar_files.txt src/data/download_tar_folders.sh
	src/data/download_tar_folders.sh $@

## 2. Clean the raw OTU tables and metadata files
# Note: technically these clean files should depend on the raw_tar_files,
# but because the stem of the files doesn't necessarily match, including
# raw_tar_files as a prerequisite means that make can no longer parallelize
# the cleaning steps. So if you changed the raw data for one of the OTU tables,
# you need to delete it so that make knows to re-process it.
# If my results_folders were better labeled, I could simply
# write a rule like %.otu_table.clean : clean.py raw_data/%.tar.gz

# This code cleans both the OTU and metadata files,
# and writes both *.otu_table.clean and *.metadata.clean
$(clean_otu_tables): src/data/clean_otu_and_metadata.py $(yaml_file)
	python src/data/clean_otu_and_metadata.py data/raw_otu_tables \
	$(yaml_file) $@

# Recover from the removal of $@
# i.e. if the metadata file is deleted but the OTU table still is unchanged
$(clean_metadata_files): $(clean_otu_tables)
	@if test -f $@; then :; else \
		rm -f $(subst metadata,otu_table,$@); \
		$(MAKE) $(AM_MAKEFLAGS) $(subst metadata,otu_table,$@); \
  	fi

## 3. Manual meta-analysis
# Maybe pull this from the internet? Maybe just have it echo the link to download it from? Or, like, "grab it from the paper"?
$(manual_meta_analysis):
	echo -e "You can find the manual meta analysis files in data/lit_search"

####################
##### ANALYSES #####
####################

## Define the target analysis files
# Note that the file names for these qvalues files are hard-coded in the
# respective scripts, rather than being passed in as inputs.
qvalues = data/analysis_results/q-val_all_results.mean.kruskal-wallis.case-control.txt
meta_qvalues = data/analysis_results/meta.counting.q-0.05.disease_wise.txt
overall_qvalues = data/analysis_results/meta.counting.q-0.05.2_diseases.across_all_diseases.txt

dataset_info = data/analysis_results/datasets_info.txt
dysbiosis = data/analysis_results/dysbiosis_metrics.txt

alpha_divs = data/analysis_results/alpha_diversity.txt
alpha_pvals = data/analysis_results/alpha_diversity.pvalues.txt

rf_results = data/analysis_results/rf_results.txt
rf_param_search = data/analysis_results/rf_results.parameter_search.txt

ubiquity = data/analysis_results/ubiquity_abundance_calculations.txt

analysis: qvals $(dysbiosis) alpha rf_results

# Make this separately, because it takes forever
rf_param_search: $(rf_param_search)

# Some other nice subsets of the analyses, mostly for testing
qvals: $(qvalues) $(meta_qvalues) $(overall_qvalues)
alpha: $(alpha_divs) $(alpha_pvals)
rf_results: $(rf_results)

## 1. q-values files for all genera across all studies
$(qvalues): src/analysis/get_qvalues.py $(clean_otu_tables) $(clean_metadata_files)
	python src/analysis/get_qvalues.py data/clean_tables data/analysis_results

## 2. meta-analysis results (+/- 1's) for each disease
$(meta_qvalues): src/analysis/meta_analyze.py $(qvalues)
	python src/analysis/meta_analyze.py $(qvalues) data/analysis_results 0.05 2

## 3. overall meta-analysis results are made at the same time as meta_qvalues
$(overall_qvalues): $(meta_qvalues)
	@if test -f $@; then :; else \
		rm -f $(meta_qvalues); \
		$(MAKE) $(AM_MAKEFLAGS) $(meta_qvalues); \
	fi

## 4. dysbiosis metrics
$(dysbiosis): src/analysis/dysbiosis_metrics.py $(qvalues) $(dataset_info) $(overall_qvalues) $(rf_results)
	python src/analysis/dysbiosis_metrics.py $(qvalues) $(dataset_info) \
	$(overall_qvalues) $(rf_results) $(dysbiosis)

## 5. alpha diversities
$(alpha_divs): src/analysis/alpha_diversity.py $(clean_otu_tables) $(clean_metadata_files)
	python src/analysis/alpha_diversity.py data/clean_tables \
	$(alpha_divs) $(alpha_pvals)

$(alpha_pvals): $(alpha_divs)
	@if test -f $@; then :; else \
		rm -f $(alpha_pvals); \
		$(MAKE) $(AM_MAKEFLAGS) $(alpha_pvals); \
	fi

## 6. random forest results
$(rf_results): src/analysis/classifiers.py $(clean_otu_tables) $(clean_metadata_files)
	python src/analysis/classifiers.py data/clean_tables \
	$(rf_results)

## 7. random forest parameter search
$(rf_param_search): src/analysis/classifiers_parameters.py $(clean_otu_tables) $(clean_metadata_files)
	python src/analysis/classifiers_parameters.py data/clean_tables \
	$(rf_param_search)

## 8. Ubiquity and abundance
ubiquity: $(ubiquity)
$(ubiquity): src/analysis/ubiquity_abundance.py $(clean_otu_tables) $(clean_metadata_files) $(overall_qvalues)
	python src/analysis/ubiquity_abundance.py data/clean_tables \
	$(overall_qvalues) $@

#######################
##### PHYLOT TREE #####
#######################

genera_file = data/analysis_results/genera.tmp
ncbi_file = data/analysis_results/ncbi_ids.tmp
clean_ncbi = data/analysis_results/ncbi_ids.clean.for_phyloT
phyloT_file = data/user_input/phyloT_tree.newick
final_tree_file = data/analysis_results/phyloT_tree.updated.newick

# Note: phyloT doesn't necessarily return genera in the same order for
# the same input list of genus IDs, so the arrangement of genera may
# not be exactly as in the paper figure. This is okay - equivalent
# trees have many different linear representations
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
	src/analysis/clean_ncbi.py $(ncbi_file) \
	data/analysis_results/ncbi_ids.clean.tmp $(clean_ncbi)

# Manual step: go to the phyloT website and make the tree
$(phyloT_file): $(clean_ncbi)
	read -n1 -p "Go to http://phylot.biobyte.de/ and generate tree from ${clean_ncbi}. Press any key to continue once you've added the tree file in ${phyloT_file}. "

# Manually edit tree with genera which didn't have NCBI IDs
$(final_tree_file): src/analysis/update_tree.py $(phyloT_file) $(genera_file)
	src/analysis/update_tree.py $(genera_file) $(phyloT_file) $(final_tree_file)

#############################
##### PREP FOR PLOTTING #####
#############################

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
	python src/analysis/reorder_qvalues.py $(qvalues) --qthresh 0.05  \
	$(meta_qvalues) $(overall_qvalues) $(final_tree_file)

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
	python src/analysis/logfold_effect.py data/clean_tables \
	$(qvalues_clean) $(logfold)

##################
##### TABLES #####
##################

## Tables
# Table 1 (main text) has the datasets, controls, cases, and references
table1 = final/tables/table1.dataset_info.tex
# Table 2 (supplement) has the same things, plus some info about the sequencers
table2 = final/tables/table2.dataset_info_supplement.tex
# Table 3 has the processing parameters
table3 = final/tables/table3.processing_info.tex
# Table 4 has the data and metadata sources
table4 = final/tables/table4.data_metadata_sources.tex

tables: $(table1) $(table2) $(table3) $(table4)

# tables-1-2.dataset_info.py makes table1, table2, and dataset_info
$(table1): src/final/table.datasets_info.py $(yaml_file) $(clean_otu_tables) $(clean_metadata_files)
	python $< $(yaml_file) data/raw_otu_tables data/clean_tables \
	$(dataset_info) $(table1) $(table2)

$(table2): $(table1)
	@if test -f $@; then :; else \
		rm -f $(table1); \
		$(MAKE) $(AM_MAKEFLAGS) $(table1); \
	fi

$(dataset_info): $(table1)
	@if test -f $@; then :; else \
		rm -f $(table1); \
		$(MAKE) $(AM_MAKEFLAGS) $(table1); \
	fi

$(table3): src/final/table.processing_info.py $(yaml_file)
	python $< $(yaml_file) data/raw_otu_tables $(table3) $(table4)

$(table4): $(table3)
	@if test -f $@; then :; else \
		rm -f $(table3); \
		$(MAKE) $(AM_MAKEFLAGS) $(table3); \
	fi
###################
##### FIGURES #####
###################

figures: main_figures supp_figures

# Some subset of figures
main_figures: figure1 figure2 figure3
supp_figures: figure4 figure5 figure6 figure7 figure8 figure9
rf_param_figures: figure10 figure11

## Define figure file names
figure1 = final/figures/figure1.samplesize_auc_extent_direction.png

# Disease-specific heatmaps
figure2 = final/figures/figure2.cdi_heatmap.png \
          final/figures/figure2.ob_heatmap.png \
		  final/figures/figure2.ibd_heatmap.png \
		  final/figures/figure2.hiv_heatmap.png \
		  final/figures/figure2.crc_heatmap.png
figure6 = final/figures/figure6.cdi_heatmap.with_labels.png \
          final/figures/figure6.ob_heatmap.with_labels.png \
		  final/figures/figure6.ibd_heatmap.with_labels.png \
		  final/figures/figure6.hiv_heatmap.with_labels.png \
		  final/figures/figure6.crc_heatmap.with_labels.png

# Core response
figure3a = final/figures/figure3a.core_disease_with_phylo.png
figure3b = final/figures/figure3b.core_overlap.png
figure3c = final/figures/figure3c.abundance.png \
           final/figures/figure3c.ubiquity.png
figure7 = final/figures/figure7.core_disease_with_phylo.with_labels.png

# Alpha diversity
figure4 = final/figures/figure4.alpha_diversity.png

# RF supplementary figures
figure5 = final/figures/figure5.roc_curves.png

# Figure 8 and 9 are the big ol' heatmaps. Not done yet.
figure8 = final/figures/figure8.overall_heatmap_log10qvalues.png
figure9 = final/figures/figure9.overall_heatmap_log2effect.png

# Note: figures 11 and 12 should NOT be in make 'all',
# they should be with rf_params
figure10 = final/figures/figure10.rf_params_gini.png
figure11 = final/figures/figure11.rf_params_entropy.png

figure1: $(figure1)
figure2: $(figure2)
figure3: $(figure3a) $(figure3b) $(figure3c)
figure4: $(figure4)
figure5: $(figure5)
figure6: $(figure6)
figure7: $(figure7)
figure8: $(figure8)
figure9: $(figure9)
figure10: $(figure10)
figure11: $(figure11)

# Figure 1: sample size, AUC, extent and direction of shifts
$(figure1): src/final/figure.samplesize_auc_extent_direction.py $(dysbiosis) $(dataset_info)
	python $< $(dysbiosis) $(dataset_info) $(figure1)

# Figure 2: disease heatmaps without labels
# $* contains the disease string, $@ is the target file
final/figures/figure2.%_heatmap.png: src/final/figure.disease_specific_heatmaps.py $(qvalues) $(dataset_info)
	python $< $* $(qvalues) $(dataset_info) $@

# Figure 3A: Core heatmaps: disease-wise, core, and phylogeny
$(figure3a): src/final/figure.core_and_disease_specific_genera.py $(meta_clean) $(overall_clean)
	python $< $(meta_clean) $(overall_clean) $@

# Figure 3B: Percent overlap
$(figure3b): src/final/figure.percent_overlap.py $(dysbiosis) $(dataset_info)
	python $< $(dysbiosis) $(dataset_info) $@

# Figure 3C: Ubiquity and abundance
final/figures/figure3c.%.png: src/final/figure.ubiquity_abundance_boxplots.py $(ubiquity)
	python $< $(ubiquity) $* $@

# Figure 4: alpha diversities
$(figure4): src/final/figure.alpha_diversity.py $(alpha_divs)
	python $< $(alpha_divs) $@

# Figure 5: ROC curves
$(figure5): src/final/figure.roc_curves.py $(rf_results)
	python $< $(rf_results) $@

# Figure 6: Disease-specific heatmap with labels
final/figures/figure6.%_heatmap.with_labels.png: src/final/figure.disease_specific_heatmaps.py $(qvalues) $(dataset_info)
	python $< $* $(qvalues) $(dataset_info) $@ --labels

# Figure 7: Core heatmap with labels
$(figure7): src/final/figure.core_and_disease_specific_genera.py $(meta_clean) $(overall_clean)
	python $< $(meta_clean) $(overall_clean) $@ --labels

# Figure 8: Overall heatmap with q values
$(figure8): src/final/figure.overall_heatmap.py $(qvalues_clean) $(meta_clean) $(overall_clean) $(dataset_info)
	python $< $(qvalues_clean) $(meta_clean) $(overall_clean) $(dataset_info) $@ --plot-log10

# Figure 9: Overall hetamap with effect
$(figure9): src/final/figure.overall_heatmap.py $(logfold) $(meta_clean) $(overall_clean) $(dataset_info)
	python $< $(logfold) $(meta_clean) $(overall_clean) $(dataset_info) $@

# Figure 10: RF parameter search, gini criteria
$(figure10): src/final/figure.rf_params.py $(rf_param_search)
	python $< $(rf_param_search) gini $@

# Figure 11: RF parameter search, entropy criteria
$(figure11): src/final/figure.rf_params.py $(rf_param_search)
	python $< $(rf_param_search) entropy $@

###############################
##### SUPPLEMENTARY FILES #####
###############################

supp_qvals = final/supp-files/file-S1.qvalues.txt
supp_disease = final/supp-files/file-S2.disease_specific_genera.txt
supp_overall = final/supp-files/file-S3.core_genera.txt
supp_litsearch = final/supp-files/file-S4.literature_results.txt
supp_files: $(supp_qvals) $(supp_disease) $(supp_overall) $(supp_litsearch)

$(supp_qvals): $(qvalues)
	cp $(qvalues) $@

$(supp_disease): src/final/supp-file.convert_meta_analysis_results.py $(meta_qvalues)
	python $< $(meta_qvalues) $@

$(supp_overall): src/final/supp-file.convert_meta_analysis_results.py $(overall_qvalues)
	python $< $(overall_qvalues) $@

$(supp_litsearch): $(manual_meta_analysis)
	cp $< $@
