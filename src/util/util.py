#!/usr/bin/env python

"""
Useful functions to be used through data processing and analysis code.
"""

import copy
import pandas as pd

# Stats functions
from scipy.stats.mstats import kruskalwallis
from scipy.stats import ranksums, mannwhitneyu
# FDR correction
from statsmodels.sandbox.stats.multicomp import multipletests

def raw2abun(df):
    """
    Converts OTU table with counts to relative abundances.

    Parameters
    ----------
    df : pandas dataframe
        OTUs in columns, samples in rows

    Returns
    -------
    df : pandas dataframe
        input dataframe normalized by total number of reads per sample
    """
    return df.divide(df.sum(axis=1), axis=0)

def collapse_taxonomic_contents_df(OTU_table, taxonomic_level):
    """
    Collapses OTU table to given taxonomic level by string-matching.

    Adapated from Thomas Gurry's code at
    https://github.com/thomasgurry/amplicon_sequencing_pipeline

    Parameters
    ----------
    OTU_table : pandas dataframe
        OTUs in columns, samples in rows.
        Taxonomic levels in OTU strings should be semicolon-delimited,
        starting with kingdom level.
        Unannotated taxonomic levels should end with '__' (e.g. ''...;g__Roseburia;s__')
    taxonomic_level : str
        kingdom, phylum, class, order, family, genus, or species

    Returns
    -------
    newdf : pandas dataframe
        OTUs in columns, samples in rows.
        OTUs are collapsed to the given taxonomic level.
        Matching values (for annotated taxa) are summed for each sample.
        Values corresponding to unannotated taxa are discarded.
    """

    OTU_IDs = list(OTU_table.columns)

    # Collapse to the right level
    if(taxonomic_level == "kingdom"):
        OTU_taxa = [OTU_ID.split(';')[0] for OTU_ID in OTU_IDs]
    if(taxonomic_level == "phylum"):
        OTU_taxa = [';'.join(OTU_ID.split(';')[:2]) for OTU_ID in OTU_IDs]
    if(taxonomic_level == "class"):
        OTU_taxa = [';'.join(OTU_ID.split(';')[:3]) for OTU_ID in OTU_IDs]
    if(taxonomic_level == "order"):
        OTU_taxa = [';'.join(OTU_ID.split(';')[:4]) for OTU_ID in OTU_IDs]
    if(taxonomic_level == "family"):
        OTU_taxa = [';'.join(OTU_ID.split(';')[:5]) for OTU_ID in OTU_IDs]
    if(taxonomic_level == "genus"):
        OTU_taxa = [';'.join(OTU_ID.split(';')[:6]) for OTU_ID in OTU_IDs]
    if(taxonomic_level == "species"):
        OTU_taxa = [';'.join(OTU_ID.split(';')[:7]) for OTU_ID in OTU_IDs]

    # Get indices of each unique taxon
    taxa_indices = {}
    for i in range(len(OTU_taxa)):
        if (OTU_taxa[i] not in taxa_indices) and (OTU_taxa[i][(len(OTU_taxa[i])-2):] != "__"):
            taxa_indices[OTU_taxa[i]] = []
        if (OTU_taxa[i][(len(OTU_taxa[i])-2):] != "__"):
            taxa_indices[OTU_taxa[i]].append(i)
    # Make new empty df with the same samples as original df and taxa in taxa_indices.keys()
    newdf = pd.DataFrame(index=OTU_table.index, columns=taxa_indices.keys(), data=0)

    # Get sample contents for each taxa of the chosen level and put into newdf
    for key in taxa_indices:
        indices = taxa_indices[key]
        newcol = OTU_table.iloc[:, indices].sum(axis=1)
        newdf[key] = copy.copy(newcol)

    return newdf


def compare_otus_teststat(df, Xsmpls, Ysmpls, method='kruskal-wallis', multi_comp=None):
    """
    Compares columns between Xsmpls and Ysmpls, with statistical method=method.
    Returns dataframe with both the qvals ('p') and test statistic ('test-stat')

    parameters
    ----------
    df             dataframe, samples are in rows and OTUs in columns
    X,Ysmpls       list of samples to compare
    method         statistical method to use for comparison
    multi_comp     str, type of multiple comparison test to do.
                   Currently accepts 'fdr' or None

    outputs
    -------
    results        dataframe with OTUs in rows and 'p' and 'test-stat' in columns

    """
    if method == 'kruskal-wallis':
        pfun = kruskalwallis
    elif method == 'wilcoxon' or method == 'ranksums':
        pfun = ranksums
    elif method == 'mann-whitney':
        pfun = mannwhitneyu
        # Note: prob wanna add some kwargs here to say whether 2sided or not

    results = pd.DataFrame(index=df.columns, columns=['test-stat', 'p'])
    for o in df.columns:
        try:
            h, p = pfun(df.loc[Xsmpls, o], df.loc[Ysmpls, o])
        except:
            p = 1
            h = 0
        results.loc[o, 'p'] = p
        results.loc[o, 'test-stat'] = h

    if multi_comp == 'fdr':
        _, results['q'], _, _ = multipletests(results['p'], method='fdr_bh')

    return results
