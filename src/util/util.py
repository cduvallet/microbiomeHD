#!/usr/bin/env python

"""
Useful functions to be used through data processing and analysis code.
"""

import copy
import pandas as pd
import numpy as np

# Stats functions
from scipy.stats.mstats import kruskalwallis
from scipy.stats import ranksums, mannwhitneyu
# FDR correction
from statsmodels.sandbox.stats.multicomp import multipletests
# Classifiers
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import auc, roc_curve, confusion_matrix, cohen_kappa_score, make_scorer
from sklearn.cross_validation import StratifiedKFold
from scipy import interp
from scipy.stats import fisher_exact

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

def prep_classifier(df, H_smpls, dis_smpls, random_state):
    """
    Prepares a classifier given a dataframe and list of samples for each
    class.

    Parameters
    ----------
    df : pandas dataframe
        dataframe with samples in rows, features in columns
    H_smpls, dis_smpls : list (or pandas index)
        samples in each class, should be in index of df
    random_state : int
        classifier seed

    Returns
    -------
    rf : RandomForestClassifier
        classifier with 1000 estimators, oob_score=True, and
        random_state seed
    X : numpy array
        values from df to be used for classification
    Y : list
        vector of 0's and 1's. 1 corresponds to samples in dis_smpls,
        0 is samples in H_smpls. Samples in X and Y are in the same
        order.
    """
    if not isinstance(H_smpls, list):
        H_smpls = list(H_smpls)
    if not isinstance(dis_smpls, list):
        dis_smpls = list(dis_smpls)

    all_smpls = H_smpls + dis_smpls

    rf = RandomForestClassifier(n_estimators=1000, random_state=random_state)
    X = df.loc[all_smpls].values
    Y = [1 if i in dis_smpls else 0 for i in all_smpls]
    return rf, X, Y

def cv_and_roc(rf, X, Y, num_cv=5, random_state=None):
    """
    Perform cross validated training and testing and return the aggregate
    interpolated ROC curve and confusion matrices.

    Parameters
    ----------
    rf : any sklearn classifier object
    X : array-like or sparse matrix, shape = [n_samples, n_features]
        The input samples to be split into train and test folds and
        cross-validated.
    Y : list or array
        array-like, shape = [n_samples] or [n_samples, n_outputs]
        The target values (class labels in classification).
    num_cv : int (default: 5)
        number of cross-validation folds
    random_state : int (default 12345)
        random state seed for StratifiedKFold

    Returns
    -------
    d : dict with the following keys:
        'roc_auc': area under the ROC curve
        'conf_mat': confusion matrix array, looks like:
                              pred
                             0   1
                    true  0  -   -
                          1  -   -
        'fisher_p': fisher exact pvalue of conf_mat above
        'y_probs': probability of being class 1
        'y_trues': true labels
        'mean_fpr', 'mean_tpr': interpolated values used to build ROC curve

    """
    if isinstance(Y, list):
        Y = np.asarray(Y)
    cv = StratifiedKFold(Y, num_cv, shuffle=True, random_state=random_state)
    mean_tpr = 0.0
    mean_fpr = np.linspace(0, 1, 100)
    conf_mat = np.asarray([[0,0],[0,0]])
    y_probs = np.empty_like(Y, dtype=float)
    y_trues = np.empty_like(Y)
    y_preds = np.empty_like(Y)
    cv_count = 0
    cv_counts = np.empty_like(Y)

    for train_index, test_index in cv:
        X_train, X_test = X[train_index], X[test_index]
        Y_train, Y_test = Y[train_index], Y[test_index]
        probs = rf.fit(X_train, Y_train).predict_proba(X_test)[:,1]

        # Store probability and true Y for later
        y_probs[test_index] = probs
        y_trues[test_index] = Y_test # literally redundant, but keep it to maintain backward compatibility
        y_pred = rf.predict(X_test)
        y_preds[test_index] = y_pred
        # Compute ROC curve and area under the curve
        fpr, tpr, thresholds = roc_curve(Y_test, probs)
        mean_tpr += interp(mean_fpr, fpr, tpr)
        # Compute confusion matrix
        conf_mat += confusion_matrix(Y_test, y_pred, labels=[0,1])

        # Track which fold each sample was tested in
        cv_counts[test_index] = cv_count
        cv_count += 1

    mean_tpr /= len(cv)
    roc_auc = auc(mean_fpr, mean_tpr)

    _, fisher_p = fisher_exact(conf_mat)

    return {i: j for i, j in
            zip(('roc_auc', 'conf_mat', 'mean_fpr', 'mean_tpr',
                'fisher_p', 'y_prob', 'y_true', 'test_fold',
                'y_preds'),
               (roc_auc, conf_mat, mean_fpr, mean_tpr, fisher_p,
               y_probs, y_trues, cv_counts, y_preds))}

def shuffle_col(col):
    """
    Shuffles the index labels of one column after dropping NaN's.

    Calling df.apply(shuffle_col) shuffles each column values but preserves
    the location of NaN's.

    (Why does this work? Because this function returns a Series with only
    the non-NaN indices. When you call shuffle_col on an entire dataframe,
    I think pandas just fills in the missing indices with NaN's...)
    """
    newcol = pd.Series(
        data=col.dropna().values,
        index=np.random.permutation(col.dropna().index))
    return newcol
