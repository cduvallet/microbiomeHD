# microbiome-meta-analysis
Cross-disease comparison of case-control gut microbiome studies

# Eventual directory structure

This structure follows what's recommended by [Cookie Cutter Data Science](https://drivendata.github.io/cookiecutter-data-science/)
Here's to hoping it works!

```
|-- License? Citation?
|-- Makefile   ## <- Makefile with commands like `make data` or `make figures` (or `make paper`?)
|-- README.md
|-- data
|    |-- raw OTU table results and metadata files
|    |-- cleaned up OTU tables
|    |-- lit search results (the manual stuff I did)
|    |-- analysis results files
|         |-- q values
|     	  |-- disease-wise meta-analysis results
|     	  |-- overall meta-analysis results
|     	  |-- random forest AUCs, fisher P values, sample size, 
|         |   prevalence (input to supp figures)
|
|-- src
|    |-- __init__.py  ## If I have an init.py will I be able to import
|    |                ## modules from here?!
|    |
|    |-- data
|    |    |-- maybe a script to download the raw OTU tables from wherever
|    |    |   they are?
|    |    |-- script to clean up the raw data into clean OTU tables and metadata
|    |	  |-- script to get info about the data
|    |
|    |-- analysis
|    |	  |-- script to make the qvalues in each dataset for each genus
|    |	  |-- random forests - AUCs and corr(AUC, sample size or prevalence)
|    |    |-- disease-wise meta-analysis results
|    |	  |-- overall meta-analysis results (unless that's easier in the
|    |    |   same script as above)
|    |	  |-- phyloT stuff
|    |	  |-- alpha diversities
|    |	  |-- supp random forest: playing with parameters (warning:
|    |    |   this takes forever)
|    |
|    |-- figures
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
|              |-- the alpha diversity figure
|              |-- the heatmap with p-values and effect sizes for all
|              |    genera in all datasets
|              |-- the random forest playing with parameters figures
|
|-- final
|    |-- figures
|    |    |-- save all of the figures here
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
