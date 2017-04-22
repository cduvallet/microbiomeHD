"""
Make table of the log-fold change for each genus in each dataset.
"""

import numpy as np
import pandas as pd
from copy import copy

def get_log_change(col, dis_smpls, H_smpls, method='mean', logfun=np.log2):
    """
    Get the log(dis_samples/H_smpls) for one Series of relative abundances.
    This function is to be used with the OTU table,
    df.apply(lambda col: get_log_change(col, dis_smpls, H_smpls, ...))

    Parameters
    ----------
    col : pandas Series
        series with genera in index
    dis_smpls, H_smpls : lists
        lists of samples to compare
    method : str
        'mean' or 'median'
    logfun : function
        np.log10 or np.log2

    Returns
    -------
    val : float
        logfun(col.loc[dis_smpls].method()/col.loc[H_smpls].method()).
        If method(dis) == 0, returns 0. If method(H) == 0, returns np.inf
    """
    if method == 'mean':
        try:
            val = logfun(col.loc[dis_smpls].mean()/col.loc[H_smpls].mean())
        except:
            # If both numerator and denominator are 0, return 0 logfold change
            if col.loc[dis_smpls].mean() == 0:
                val = 0
            # If numerator is not zero, return log10(1/0) AKA -log10(0/1) = np.inf
            else:
                val = np.inf

    elif method == 'median':
        try:
            val = logfun(df.loc[dis_smpls, o].median()/df.loc[H_smpls, o].median())
        except:
            # If both numerator and denominator are 0, return 0 logfold change
            if df.loc[dis_smpls, o].median() == 0:
                val = 0
            # If numerator is not zero, return log10(1/0) AKA -log10(0/1) = np.inf
            else:
                val = np.inf
    else:
        raise ValueError('Unrecognized method to compare values')
    return val

def convert_dataset_to_logfold(col, dfdict, logfun=np.log2, method='mean'):
    """
    Calculate the logfold change between disease and controls for genera
    in one column.

    Parameters
    ----------
    col : pandas Series
        genera in rows, dataset as Series name. Values are irrelevant.
    dfdict : dict
        {dataset: {'df': df, 'dis_smpls': dis_smpls, 'H_smpls': H_smpls}, ...}
    logfun : function
        np.log10 or np.log2
    method : str
        'mean' or 'median'

    Returns
    -------
    logcol : pandas Series
        Same index as col. For genera which are present in that dataset's
        OTU table, the log-fold change is filled in. For genera which are
        not present in that dataset's OTU table, value is np.nan.
        Genera which are in the OTU table but not in col.index are discarded.
    """
    dataset = col.name
    df = dfdict[dataset]['df']
    dis_smpls = dfdict[dataset]['dis_smpls']
    h_smpls = dfdict[dataset]['H_smpls']
    # This returns a Series with genera as rows, logfold change as values
    logdf = df.apply(lambda col: get_log_change(col, dis_smpls, h_smpls,
                                                logfun=logfun, method=method))

    # Merge results from logfold change with the given index, and keep only
    # rows which were in the original index. Any rows which are not in the df
    # are nan's
    logcol = pd.Series(index=col.index, data=np.nan)

    genera = [i for i in logcol.index if i in logdf.index]

    logcol.loc[genera] = logdf.loc[genera]

    return logcol


# Set up arguments, including logfunction and method (mean or median)
# and input files: datadir, qvalues to be converted (the reordered ones prob)

# Read in dfdict

# Collapse to genus level

# Read in qvalues

# Calculate logfold
allres = qvals.apply(lambda col: convert_dataset_to_logfold(col, dfdict, logfun=np.log2, method='mean'))

# Replace +/- infinity with max/min value in entire matrix
# From get_log_change():
# +infinity is returned to controls = 0, disease > 0
# -infinity is returned when controls > 0, disease = 0
# (note: 0 is returned when both controls and disease = 0)
allres = allres.replace(np.inf, np.ma.masked_invalid(allres.fillna(0)).max())
allres = allres.replace(-np.inf, np.ma.masked_invalid(allres.fillna(0)).min())

# Write results
