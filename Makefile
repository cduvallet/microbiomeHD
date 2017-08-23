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
reviewer_comments: reviewer_analysis

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

dataset_info = data/analysis_results/datasets_info.txt
split_dataset_info = data/analysis_results/datasets_info.split_cases.txt

raw_data: $(raw_tar_files)
clean_data: $(clean_otu_tables) $(clean_metadata_files)
data: raw_data clean_data $(manual_meta_analysis) $(dataset_info)

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

## 4. Dataset info - table with basic information about the datasets
$(dataset_info): src/data/dataset_info.py $(yaml_file) $(clean_otu_tables) $(clean_metadata_files)
	python $< $(yaml_file) data/raw_otu_tables data/clean_tables $@

# Same as above, but with case patients split into separate groups
$(split_dataset_info): src/data/dataset_info.py $(yaml_file) $(clean_otu_tables) $(clean_metadata_files) $(split_datasets)
	python $< $(yaml_file) data/raw_otu_tables data/clean_tables $@ --split-cases --subset $(split_datasets)

####################
##### ANALYSES #####
####################

## Define the target analysis files
# Note that the file names for these qvalues files are hard-coded in the
# respective scripts, rather than being passed in as inputs.
qvalues = data/analysis_results/qvalues.mean.kruskal-wallis.case-control.txt
meta_qvalues = data/analysis_results/meta.counting.q-0.05.disease_wise.txt
# Based on counting heuristic (significant in at least 2 diseases)
overall_qvalues = data/analysis_results/meta.counting.q-0.05.2_diseases.across_all_diseases.txt
# Based on counting heuristic excluding diarrhea datasets
nocdi_overall = data/analysis_results/meta.counting.q-0.05.2_diseases.across_all_diseases_except_cdi.txt
# +/- 1 labels, based on combining qvalues with Stouffer's method
overall_qvalues_stouffer = data/analysis_results/meta.stouffer.q-0.05.across_all_diseases.txt
# Raw combined pvalues from Stouffer's method
qvalues_stouffer = data/analysis_results/meta.stouffer_qvalues.txt

dysbiosis = data/analysis_results/dysbiosis_metrics.txt

alpha_divs = data/analysis_results/alpha_diversity.txt
alpha_pvals = data/analysis_results/alpha_diversity.pvalues.txt

rf_results = data/analysis_results/rf_results.txt
rf_core = data/analysis_results/rf_results.core_only.txt
rf_param_search = data/analysis_results/rf_results.parameter_search.txt
rf_h_v_dis = data/analysis_results/rf_results.healthy_vs_disease.txt

ubiquity = data/analysis_results/ubiquity_abundance_calculations.txt

# Reviewer comment: split analysis for UC/CD case patients
split_datasets = data/user_input/split_cases_datasets.txt
split_qvalues = data/analysis_results/qvalues.mean.kruskal-wallis.split-cases.txt
split_rf = data/analysis_results/rf_results.split_cases.txt
split_dysbiosis = data/analysis_results/dysbiosis_metrics.split_cases.txt

analysis: qvals $(dysbiosis) alpha rf_results
reviewer_analysis: $(rf_core) $(overall_qvalues_stouffer) $(stouffer) $(nocdi_overall) $(split_qvalues) $(split_dysbiosis) $(split_rf)

# Make this separately, because it takes forever
rf_param_search: $(rf_param_search)

# Some other nice subsets of the analyses, mostly for testing
qvals: $(qvalues) $(meta_qvalues) $(overall_qvalues) $(split_qvalues)
alpha: $(alpha_divs) $(alpha_pvals)
rf_results: $(rf_results) $(rf_h_v_dis)

## 1. q-values files for all genera across all studies
$(qvalues): src/analysis/get_qvalues.py $(clean_otu_tables) $(clean_metadata_files)
	python $< data/clean_tables $@

## 2. meta-analysis results (+/- 1's) for each disease
$(meta_qvalues): src/analysis/meta_analyze.py $(qvalues)
	python $< $(qvalues) data/analysis_results 0.05 2 --exclude-nonhealthy

## 3. overall meta-analysis results are made at the same time as meta_qvalues
$(overall_qvalues): $(meta_qvalues)
	@if test -f $@; then :; else \
		rm -f $(meta_qvalues); \
		$(MAKE) $(AM_MAKEFLAGS) $(meta_qvalues); \
	fi

# Overall meta-analysis with heuristic, excluding diarrhea datasets
$(nocdi_overall): src/analysis/meta_analyze.py $(qvalues)
	python $< --no-cdi $(qvalues) data/analysis_results 0.05 2

# Overall meta-analysis using Stouffer's method for combining pvalues
$(overall_qvalues_stouffer): src/analysis/meta_analyze_stouffer.py $(qvalues) $(dataset_info)
	python $< $(qvalues) $(dataset_info) $(qvalues_stouffer) $@

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

## 9. Random forest using only core bugs
$(rf_core): src/analysis/classifiers.py $(clean_otu_tables) $(clean_metadata_files) $(overall_qvalues)
	python $< --core $(overall_qvalues) data/clean_tables $@

## 10. Random forest for general healthy vs disease classifier
$(rf_h_v_dis): src/analysis/healthy_disease_classifier.py $(clean_otu_tables) $(clean_metadata_files)
	python $< data/clean_tables $@

## 10. Reviewer comment: re-do major analyses for subgroups of case patients
## separately
$(split_qvalues): src/analysis/get_qvalues.py $(split_datasets) $(clean_otu_tables) $(clean_metadata_files)
	python $< data/clean_tables $@ --subset $(split_datasets) --split-cases

$(split_rf): src/analysis/classifiers.py $(split_datasets) $(clean_otu_tables) $(clean_metadata_files)
	python $< data/clean_tables $@ --subset $(split_datasets) --split-cases

$(split_dysbiosis): src/analysis/dysbiosis_metrics.py $(split_qvalues) $(split_dataset_info) $(overall_qvalues) $(split_rf)
	python $< $(split_qvalues) $(split_dataset_info) \
	$(overall_qvalues) $(split_rf) $@

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
qvalues_clean_tmp = $(subst txt,sig.txt,$(qvalues))
qvalues_clean = $(subst txt,sig_ordered.txt,$(qvalues))
meta_clean = $(subst txt,sig_ordered.txt,$(meta_qvalues))
overall_clean = $(subst txt,sig_ordered.txt,$(overall_qvalues))
logfold = $(subst txt,log2change.sig_ordered.txt,$(qvalues))

# Different core definitions
nocdi_clean = $(subst txt,sig_ordered.txt,$(nocdi_overall))
stouffer_clean = $(subst txt,sig_ordered.txt,$(overall_qvalues_stouffer))

for_plotting: $(qvalues_clean) $(meta_clean) $(overall_clean) $(logfold)

$(qvalues_clean_tmp): src/analysis/clean_qvalues.py $(qvalues)
	python $< $(qvalues)

$(qvalues_clean): src/analysis/reorder_qvalues.py $(qvalues_clean_tmp) $(final_tree_file)
	python $< --do-qvals --qvalues $(qvalues_clean_tmp) $(final_tree_file)

$(meta_clean): src/analysis/reorder_qvalues.py $(meta_qvalues) $(final_tree_file) $(qvalues_clean)
	python $< --disease-df $(meta_qvalues) --qvalues $(qvalues_clean) $(final_tree_file)

$(overall_clean): src/analysis/reorder_qvalues.py $(overall_qvalues) $(final_tree_file) $(qvalues_clean)
	python $< --overall $(overall_qvalues) --qvalues $(qvalues_clean) $(final_tree_file)

$(nocdi_clean): src/analysis/reorder_qvalues.py $(nocdi_overall) $(final_tree_file)
	python $< --overall $(nocdi_overall) $(final_tree_file)

$(stouffer_clean): src/analysis/reorder_qvalues.py $(overall_qvalues_stouffer) $(final_tree_file)
	python $< --overall $(overall_qvalues_stouffer) $(final_tree_file)

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
$(table1): src/final/table.datasets_info.py $(dataset_info)
	python $< $(dataset_info) $(table1) $(table2)

$(table2): $(table1)
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
supp_figures: figure4 figure5 figure6 figure7 figure8 figure9 figure12 figure13 figure14 figure15
rf_param_figures: figure16 figure17

## Define figure file names
figure1 = final/figures/figure1.samplesize_auc_extent_direction.pdf
figure7 = final/figures/figure7.split_cases.samplesize_auc_extent_direction.pdf

# Disease-specific heatmaps
figure2 = final/figures/figure2.cdi_heatmap.pdf \
          final/figures/figure2.ob_heatmap.pdf \
		  final/figures/figure2.ibd_heatmap.pdf \
		  final/figures/figure2.hiv_heatmap.pdf \
		  final/figures/figure2.crc_heatmap.pdf
figure5 = final/figures/figure5.cdi_heatmap.with_labels.pdf \
          final/figures/figure5.ob_heatmap.with_labels.pdf \
		  final/figures/figure5.ibd_heatmap.with_labels.pdf \
		  final/figures/figure5.hiv_heatmap.with_labels.pdf \
		  final/figures/figure5.crc_heatmap.with_labels.pdf
figure8 = final/figures/figure8.cd_heatmap.with_labels.pdf \
          final/figures/figure8.uc_heatmap.with_labels.pdf
# Core response
figure3a = final/figures/figure3a.core_disease_with_phylo.pdf
figure3b = final/figures/figure3b.core_overlap.pdf
figure3c = final/figures/figure3c.abundance.pdf \
           final/figures/figure3c.ubiquity.pdf
figure6 = final/figures/figure6.core_disease_with_phylo.with_labels.pdf

# Alpha diversity
figure9 = final/figures/figure9.alpha_diversity.shannon.pdf

# RF supplementary figures
figure4 = final/figures/figure4.roc_curves.pdf

# The big ol' heatmaps.
figure14 = final/figures/figure14.overall_heatmap_log10qvalues.pdf
figure15 = final/figures/figure15.overall_heatmap_log2effect.pdf

# Note: figures 16 and 17 should NOT be in make 'all',
# they should be with rf_params
figure16 = final/figures/figure16.rf_params_gini.pdf
figure17 = final/figures/figure17.rf_params_entropy.pdf

# Need to figure out where these fit in...
rf_dataset_out = final/figures/figure13.rf_healthy_disease.dataset_out.pdf
rf_disease_out = final/figures/figure13.rf_healthy_disease.disease_out.pdf
healthy_dis_rf: $(rf_dataset_out) $(rf_disease_out)

# Different ways to define core bugs
core_defns_fig = final/figures/figure12.different_core_definitions.pdf
core_defns: $(core_defns_fig)

figure1: $(figure1)
figure2: $(figure2)
figure3: $(figure3a) $(figure3b) $(figure3c)
figure4: $(figure4)
figure5: $(figure5)
figure6: $(figure6)
figure7: $(figure7)
figure8: $(figure8)
figure9: $(figure9)
# figures 10 and 11 are the other alpha diversity ones
figure12: core_defns
figure13: healthy_dis_rf
figure14: $(figure14)
figure15: $(figure15)
figure16: $(figure16)
figure17: $(figure17)

## Figure dependencies
fmt = src/util/Formatting.py

# Sample size, AUC, extent and direction of shifts
$(figure1): src/final/figure.samplesize_auc_extent_direction.py $(dysbiosis) $(dataset_info) $(fmt)
	python $< $(dysbiosis) $(dataset_info) $(figure1)

$(figure7): src/final/figure.samplesize_auc_extent_direction.py $(split_dysbiosis) $(split_dataset_info) $(fmt)
	python $< $(split_dysbiosis) $(split_dataset_info) $@ --edd

# Disease heatmaps without labels
# $* contains the disease string, $@ is the target file
final/figures/figure2.%_heatmap.pdf: src/final/figure.disease_specific_heatmaps.py $(qvalues) $(dataset_info) $(fmt)
	python $< $* $(qvalues) $(dataset_info) $@

# Core heatmaps: disease-wise, core, and phylogeny
$(figure3a): src/final/figure.core_and_disease_specific_genera.py $(meta_clean) $(overall_clean) $(fmt)
	python $< $(meta_clean) $(overall_clean) $@

# Percent overlap
$(figure3b): src/final/figure.percent_overlap.py $(dysbiosis) $(dataset_info) $(fmt)
	python $< $(dysbiosis) $(dataset_info) $@

# Ubiquity and abundance
final/figures/figure3c.%.pdf: src/final/figure.ubiquity_abundance_boxplots.py $(ubiquity)
	python $< $(ubiquity) $* $@

# Alpha diversities
$(figure9): src/final/figure.alpha_diversity.py $(alpha_divs) $(fmt)
	python $< $(alpha_divs) $(subst .shannon.pdf,,$@)

# ROC curves
$(figure4): src/final/figure.roc_curves.py $(rf_results) $(fmt)
	python $< $(rf_results) $@

# Disease-specific heatmap with labels
final/figures/figure5.%_heatmap.with_labels.pdf: src/final/figure.disease_specific_heatmaps.py $(qvalues) $(dataset_info) $(fmt)
	python $< $* $(qvalues) $(dataset_info) $@ --labels

final/figures/figure8.%_heatmap.with_labels.pdf: src/final/figure.disease_specific_heatmaps.py $(split_qvalues) $(split_dataset_info) $(fmt)
	python $< $* $(split_qvalues) $(split_dataset_info) $@ --labels

# Core heatmap with labels
$(figure6): src/final/figure.core_and_disease_specific_genera.py $(meta_clean) $(overall_clean) $(fmt)
	python $< $(meta_clean) $(overall_clean) $@ --labels

# Overall heatmap with q values
$(figure14): src/final/figure.overall_heatmap.py $(qvalues_clean) $(meta_clean) $(overall_clean) $(dataset_info) $(fmt)
	python $< $(qvalues_clean) $(meta_clean) $(overall_clean) $(dataset_info) $@ --plot-log10

# Overall heatmap with effect
$(figure15): src/final/figure.overall_heatmap.py $(logfold) $(meta_clean) $(overall_clean) $(dataset_info) $(fmt)
	python $< $(logfold) $(meta_clean) $(overall_clean) $(dataset_info) $@

# RF parameter search, gini criteria
$(figure16): src/final/figure.rf_params.py $(rf_param_search) $(fmt)
	python $< $(rf_param_search) gini $@

# RF parameter search, entropy criteria
$(figure17): src/final/figure.rf_params.py $(rf_param_search) $(fmt)
	python $< $(rf_param_search) entropy $@

# Healthy vs disease classifier
$(rf_dataset_out): src/final/figure.healthy_vs_disease_classifier.py $(rf_results) $(rf_h_v_dis)
	python $< $(rf_results) $(rf_h_v_dis) $(rf_dataset_out) $(rf_disease_out)

$(rf_disease_out): $(rf_dataset_out)
	@if test -f $@; then :; else \
		rm -f $<; \
		$(MAKE) $(AM_MAKEFLAGS) $<; \
	fi

# Different core definitions
$(core_defns_fig): src/final/figure.core_different_definitions.py $(overall_clean) $(nocdi_clean) $(stouffer_clean) $(final_tree_file)
	python $< --labels $(overall_clean) $(nocdi_clean) $(stouffer_clean) $(final_tree_file) $@

###############################
##### SUPPLEMENTARY FILES #####
###############################

supp_qvals = final/supp-files/file-S1.qvalues.txt
supp_disease = final/supp-files/file-S2.disease_specific_genera.txt
supp_overall = final/supp-files/file-S3.core_genera.txt
supp_litsearch = final/supp-files/file-S4.literature_results.txt
supp_effects = final/supp-files/file-S5.effects.txt
supp_files: $(supp_qvals) $(supp_disease) $(supp_overall) $(supp_litsearch) $(supp_effects)

$(supp_qvals): $(qvalues)
	cp $(qvalues) $@

$(supp_disease): src/final/supp-file.convert_meta_analysis_results.py $(meta_qvalues)
	python $< $(meta_qvalues) $@

$(supp_overall): src/final/supp-file.convert_meta_analysis_results.py $(overall_qvalues)
	python $< $(overall_qvalues) $@

$(supp_litsearch): $(manual_meta_analysis)
	cp $< $@

$(supp_effects): $(logfold)
	cp $< $@
