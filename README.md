# microbiomeHD: Gut Microbiome in Health and Disease
A cross-disease analysis of case-control gut microbiome studies
in health and disease.

# Eventual directory structure

This structure follows what's recommended by [Cookie Cutter Data Science](https://drivendata.github.io/cookiecutter-data-science/)
Here's to hoping it works!

 ```
|-- License? Citation?
|-- Makefile   ## <- Makefile with commands like `make data` or `make figures` (or `make paper`?)
|-- README.md
|-- data
|    |-- raw_otu_tables: raw OTU table results and metadata files
|    |-- clean_tables: cleaned up OTU tables
|    |-- lit_search: lit search results (the manual stuff I did)
|    |-- user_input: user-inputted files like the results_folder.yaml and
|    |   the phyloT tree
|    |
|    |-- analysis_results
|         |-- q values
|     	  |-- disease-wise and overall meta-analysis results
|     	  |-- logfold change effects
|     	  |-- random forest AUCs, fisher P values, sample size,
|         |   prevalence (input to supp figures)
|         |-- alpha diversity results
|         |-- PhyloT tree: intermediate files (minus the user-inputted
|         |   one) and final tree
|
|-- src
|    |-- __init__.py  
|    |
|    |-- data
|    |    |-- download/copy raw results into data/raw_otu_tables
|    |    |-- script to clean up the raw data into clean OTU tables and metadata
|    |
|    |-- analysis
|    |	  |-- calculate the qvalues in each dataset for each genus
|    |    |-- do the disease-wise and overall meta-analysis analysis
|    |	  |-- random forest analysis (basic AUC, for main text figure)
|    |    |-- random forest parameter "search" (for supp. figure)
|    |	  |-- phyloT scripts to read genera, search against NCBI, and re-order
|    |    |   inputted tree with manually-added missing genera
|    |	  |-- alpha diversities
|    |	  |-- make files needed for plotting (e.g. logfold change, re-ordered
|    |    |   and only significant q-values)
|    |
|    |-- figures-tables
|         |-- main
|         |    |-- figure 1 should probably be its own thing
|         |    |-- the disease heatmaps (Fig 2 and labeled)
|         |    |-- the meta-analysis heatmap with phylogeny (Fig 3A and
|         |    |   labeled one)
|         |    |-- the core overlap figure is small but separate (Fig 3B)
|         |    |-- whatever other figure I have on here, ubiquity business
|         |        (Fig 3C)
|         |
|         |-- supplementary
|         |    |-- the alpha diversity figure
|         |    |-- the heatmap with p-values and effect sizes for all
|         |    |    genera in all datasets
|         |    |-- the random forest playing with parameters figures
|         |
|         |-- tables
|         |    |-- script to make latex tables
|-- final
|    |-- figures
|    |    |-- save all of the figures here
|    |    |-- files corresponding to each table (md and tex)
|    |
|    |-- manuscript
|         |-- README.md with link to the manuscript repo
|
|-- docs??
|    |-- how to re-run the data?
|    |-- what the Makefile does, in English? I dunno.
|        We'll see how much energy I have left after this.
|
|-- requirements.txt  ## Damnit I'm gonna have to make this aren't I? augh.
```
