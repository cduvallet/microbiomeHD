#!/usr/bin/env python
"""
This script calculates various "dysbiosis metrics" for each dataset,
given the raw data and the q-values.
"""

import argparse
import copy

import pandas as pd
import numpy as np
from scipy.stats import combine_pvalues

# Add this repo to the path
import os, sys
src_dir = os.path.normpath(os.path.join(os.getcwd(), 'src/util'))
sys.path.insert(0, src_dir)


def convert_to_one_tailed(longpvals):
    """
    Convert signed two-tailed pvalues in tidy dataframe to two columns
    of one-tailed pvalues.

    Parameters
    ----------
    longpvals : pandas dataframe
        tidy df with column 'p', containing signed p-values, where the
        sign indicates effect direction

    Returns
    -------
    longpvals : pandas dataframe
        modified in place, with additional 'p-h', and 'p-dis' columns, where
        columns contain p_value/2 (if the original p-value was in its direction)
        or 1-pvalue/2 (if the original p-value was in the other direction).
            p-h has pvalues that were negative
            p-dis has pvalues that were positive

    """
    higher_in_dis = longpvals[longpvals['p'] > 0].index
    longpvals.loc[higher_in_dis, 'p-dis'] = longpvals.loc[higher_in_dis, 'p']/2

    higher_in_h = longpvals[longpvals['p'] <= 0].index
    longpvals.loc[higher_in_h, 'p-h'] = abs(longpvals.loc[higher_in_h, 'p']/2)

    def p_for_other_side(row, side, otherside):
        if np.isnan(row[side]):
            return 1-row[otherside]
        else:
            return row[side]
    longpvals['p-dis'] = longpvals.apply(p_for_other_side,
                                         args=('p-dis', 'p-h'),
                                         axis=1)
    longpvals['p-h'] = longpvals.apply(p_for_other_side,
                                       args=('p-h', 'p-dis'),
                                       axis=1)
    return longpvals

def update_reproducibility_df_lists(values, variables, metric_labels,
                                    diseaselabels, newvalues, newvariables,
                                    newmetric_label, newdisease_label):
    """
    Update the lists that will produce the reproducibility dataframe

    Parameters
    ----------
    values, variables, metric_labels, diseaselabels : lists
        existing lists to append corresponding new values to
    newvalues : list
        list of new values to append to values,
        e.g. number of OTUs, reproducibility, balance, etc
    newvariables : list
        list of labels for values, e.g. name of OTUs, name of datasets, etc.
        Should be the same length as values.
    newmetric_label : str
        label for metric type (e.g. 'balance', 'n_sig', 'rep_score')
    newdisease_label : str
        disease label (e.g. 'cdi', 'crc', etc)

    Returns
    -------
    values, variables, metric_labels, diseaselabels : lists
        modified lists with corresponding new values appended
    """
    values += newvalues
    variables += newvariables
    metric_labels += [newmetric_label]*len(newvalues)
    diseaselabels += [newdisease_label]*len(newvalues)

    return values, variables, metric_labels, diseaselabels

def reproducibility_from_fisher(disdf, samplesizes, qthresh):
    """
    Returns the number of 'reproducible' OTUs based on weighted Fisher's method.
    Note: if I ever want to actually use these Fisher p-values, I could
    break up this function to return the `metap` dataframe

    Parameters
    ----------
    disdf : pandas dataframe
        genera in rows, datasets in columns, signed q-values in values
    samplesizes : pandas dataframe
        datasets in rows and at least column 'total' with total number
        of samples in each dataset, to use as weight for Stouffer's method
    qthresh : float
        threshold for calling a fisher meta-q value "significant"

    Returns
    -------
    n_sig : int
        total number of genera significant via fisher's method
    """

    ## Turn disdf into tidy dataframe
    longpvals = copy.deepcopy(disdf)
    longpvals['otu'] = longpvals.index
    longpvals = pd.melt(longpvals, id_vars='otu',
                        value_name='p', var_name='study')

    ## Convert two-tailed signed p-values into one-tailed pvalues
    longpvals = convert_to_one_tailed(longpvals).dropna()
    longpvals = pd.melt(longpvals, id_vars=['otu', 'study'],
                        value_vars=['p-dis', 'p-h'], var_name='pval_direction')

    ## Add sample size for each study
    longpvals['sample_size'] = \
        longpvals.apply(lambda row: samplesizes.loc[row['study'], 'total'],
                        axis=1)

    ## Get the combined p-value using weighted stouffer's method
    metap = []
    for grp, subdf in longpvals.groupby(['otu', 'pval_direction']):
        # Only consider genera which are in more than one study
        if subdf.shape[0] > 1:
            # grp is the tuple that defines the group: (otu, direction)
            direction = grp[1]
            otu = grp[0]
            numstudies = subdf.shape[0]
            # Stouffer's weight z-score test
            z, p = combine_pvalues(subdf['value'].astype(float),
                                   method='stouffer',
                                   weights=subdf['sample_size'].apply(np.sqrt))
            metap.append([otu, direction, z, p, numstudies])
    metap = pd.DataFrame(metap, columns=['otu', 'direction', 'z', 'p', 'num_studies'])

    ## Count number of significant healthy and disease bugs
    # Note that from manual inspection, it doesn't look like any genera
    # are returned as significant in both directions from this method...
    sig_h = metap.query('direction == "p-h"').query('p < @qthresh')
    sig_dis = metap.query('direction == "p-dis"').query('p < @qthresh')

    return sig_h.shape[0] + sig_dis.shape[0]

def get_dysbiosis_metrics(diseases, datasets, df, pthresh, samplesizes,
                          overall=None):
    """
    Calculate "dysbiosis" metrics in various ways.
    Metrics can be disease-wise (i.e. one number per disease) or dataset-wise
    (i.e. one number per dataset). There's also one genus-wise metric.

    Disease-wise metrics include:
        * 'rep_twostudies': the number of genera which are significant in
           the same direction in at least 2 studies within a disease.
           Note that this works by a simple sum: so if a genus is health-
           associated in 2 studies and disease-associated in 1 study, these
           cancel out and the genus is considered health-associated in 1 study
           (and thus not significant in the same direction in at least two
           studies).
        * 'rep_twostudies_norm': same number as above divided by the total
           number of studies in that disease
        * 'rep_stouffer': number of significant genera after combining
           q-values across all datasets of one disease using Stouffer's method
           for combining p-values (weighting by sample size). Note that given
           q-values are assumed to be two-tailed and are converted to
           one-tailed q-values for their given direction (indicated by the
           sign of the input q-values).
        *  'rep_stouffer_norm': same number as above divided by the total
           number of studies in that disease

    Dataset-wise metrics include:
        *  'n_sig': total number significant genera
        *  'balance': number of significant disease-associated (i.e. positive
            q-value) divided by total number of significant genera.
        *  'rep_dataset': number of genera in that study which were significant
           in the same direction in at least one other study of the same
           disease. This calculation also cancels out genera with non-matching
           directionality (i.e. a genus that is health-associated in 1 study
           and disease-associated in 2 studies is equvailent to being
           health-associated in 1 study).
        *  'rep_dataset_norm': same number as above, divided by the total
           number of significant genera in that dataset.

    The genus-wide metric, which I don't think I ever used, is:
        * 'rep_score': genera are given a score (+/- 1 if significant,
           +/- 0.5 if not significant) and summed across all datasets
           of the disease. The returned metric is this sum divided by
           the number of datasets.

    The paper (as of first submission) only actually uses 'n_sig' and 'balance'
    because these are the ones that make most sense.

    Parameters
    ----------
    diseases : list
        Diseases to consider. Should be the disease strings as in dataset IDs
    datasets : list
        list of dataset IDs to consider.
    df : pandas dataframe
        signed q-values according to effect direction. Datasets in columns,
        genera in rows. Note that 'edd_singh' should be changed to 'cdi_singh'
        if you want it considered with the other diarrhea datasets.
    pthresh : float
        significance threshold
    samplesizes : pandas dataframe
        datasets in index, ['total'] column with total sample size.
        Used in weighted fisher p-val calculation.
    overall : pandas Series
        genera in index and +/-1 in values.
        -1 means overall significant in healthy patients.
        +1 means overall significant in disease patients.

    Returns
    -------
    results : tidy pandas dataframe
        columns are: ['value', 'label', 'metric', 'disease', 'dataset'],
        where 'label' is either the disease, datasetID or the genus name
        (depending on the type of metric)
    """

    # Keep only OTUs which were signficant in at least one study
    # if x is zero, this returns zero (i.e. if there is no effect, it
    # doesn't count as significant so don't worry)
    sigmap = lambda x: np.sign(x) if abs(x) < pthresh else 0
    # This one is for the score-based metric
    simplesigmap = lambda x: np.sign(x) if abs(x) < pthresh else 0.5*np.sign(x)

    results = [[],[],[],[]]

    for dis in diseases:
        print(dis)
        ## Prepare subset df
        keep_datasets = [i for i in datasets if i.startswith(dis + '_')]
        disdf = df[keep_datasets]
        # Keep only genera which are significant in at least one study
        disdf = disdf.loc[disdf.applymap(sigmap).apply(abs).sum(axis=1) != 0]

        if disdf.empty:
            print('\tempty, everything is zero')
            metrics = ['rep_score', 'rep_twostudies', 'rep_twostudies_norm',
                       'rep_stouffer', 'rep_stouffer_norm', 'n_sig', 'balance',
                       'rep_dataset', 'rep_dataset_norm']
            for metric in metrics:
                val = 0
                if metric == "balance" or metric == 'rep_dataset':
                    val = np.nan
                if metric in ['rep_score']:
                    # No genera are significant, therefore none are scored...
                    pass
                elif metric in ['n_sig', 'balance', 'rep_dataset',
                                'rep_dataset_norm']:
                    # Dataset-wise results, each disdf column gets its own label
                    results = update_reproducibility_df_lists(
                                  *results,
                                  newvalues=[val]*disdf.shape[1],
                                  newvariables=disdf.columns.tolist(),
                                  newmetric_label=metric,
                                  newdisease_label=dis)
                elif metric in ['rep_twostudies', 'rep_twostudies_norm',
                                'rep_stouffer', 'rep_stouffer_norm']:
                    # Disease-wise results, each disease gets just one result
                    results = update_reproducibility_df_lists(
                                  *results,
                                  newvalues=[val],
                                  newvariables=[dis],
                                  newmetric_label=metric,
                                  newdisease_label=dis)

        else:
            ## Reproducibility score: +1/-1 if significant,
            # +/- 0.5 if not significant - don't weight by sample size
            # Metric is the row sum divided by number of columns
            # (i.e. sum across datasets / number of datasets), and is
            # genus-wise
            df_simple_rep = disdf.applymap(simplesigmap)
            reproducibility = list(
                df_simple_rep.sum(axis=1)/float(df_simple_rep.shape[1]))

            results = update_reproducibility_df_lists(
                          *results,
                          newvalues=reproducibility,
                          newvariables=disdf.index.tolist(),
                          newmetric_label='rep_score',
                          newdisease_label=dis)

            ## Reproducibility co-occurence: genus is 'reproducibly significant'
            # if it's sig in at least 2 studies
            # Metric returns  one number per disease)
            dfcooccur = disdf.applymap(sigmap)
            # Genus is reproducible if it is significant in the same
            # direction in at least net 2 studies
            reproducibility = sum(dfcooccur.sum(axis=1).apply(abs) > 1)
            results = update_reproducibility_df_lists(
                          *results,
                          newvalues=[reproducibility],
                          newvariables=[dis],
                          newmetric_label='rep_twostudies',
                          newdisease_label=dis)

            ## Normalize co-occurence: same as above, but value is normalized by total number of sig OTUs in that disease (one number per disease)
            reproducibility = reproducibility/float(dfcooccur.shape[0])
            results = update_reproducibility_df_lists(
                          *results,
                          newvalues=[reproducibility],
                          newvariables=[dis],
                          newmetric_label='rep_twostudies_norm',
                          newdisease_label=dis)

            ## Fisher's method: number of 'reproducible' genera, i.e. genera
            # with combine p-value < pthresh
            # This returns one number per disease
            reproducibility = \
                reproducibility_from_fisher(disdf, samplesizes, pthresh)
            results = update_reproducibility_df_lists(
                          *results,
                          newvalues=[reproducibility],
                          newvariables=[dis],
                          newmetric_label='rep_stouffer',
                          newdisease_label=dis)

            ## Normalized Fisher's method: same as above, but normalized
            # by total number of datasets in that disease
            reproducibility = reproducibility/float(disdf.shape[0])
            results = update_reproducibility_df_lists(
                          *results,
                          newvalues=[reproducibility],
                          newvariables=[dis],
                          newmetric_label='rep_stouffer_norm',
                          newdisease_label=dis)

            ## Total number of significant OTUs
            n_otus = disdf.applymap(sigmap)\
                          .applymap(lambda x: 1 if x == 1 or x == -1 else 0)\
                          .sum()
            # index of this Series is the datasets (i.e. columns of disdf)
            results = update_reproducibility_df_lists(
                          *results,
                          newvalues=list(n_otus),
                          newvariables=list(n_otus.index),
                          newmetric_label='n_sig',
                          newdisease_label=dis)

            ## Balance metric is number of significant disease-associated
            # (i.e. positive q-value) genera divided by total number of
            # significant genera
            n_pos = (disdf.applymap(sigmap) == 1).sum().astype(float)
            # If there are no significant OTUs, this needs to return nan,
            # otherwise the zero makes it look the same as "all significant
            # genera are health-associated"
            n_otus = n_otus.replace(0, np.nan)
            balance = n_pos/n_otus
            results = update_reproducibility_df_lists(
                          *results,
                          newvalues=list(balance),
                          newvariables=list(balance.index),
                          newmetric_label='balance',
                          newdisease_label=dis)

            ## Reproducibility per dataset. For each dataset, count the
            # number of significant genera which are significant (in same dir)
            # in at least one other dataset of that disease
            for dataset in disdf:
                # Keep just the genera which are significant in the one study
                coldf = disdf[dataset].apply(sigmap)
                coldf = coldf[coldf != 0]
                # If there's at least one significant bug in that study
                if coldf.shape[0] > 0:
                    coldf = disdf.loc[coldf.index].applymap(sigmap)
                    # Count how many genera are significant in at least
                    # two studies
                    reproducibility = (coldf.sum(axis=1).apply(abs) >= 2).sum()
                    results = update_reproducibility_df_lists(
                                  *results,
                                  newvalues=[reproducibility],
                                  newvariables=[dataset],
                                  newmetric_label='rep_dataset',
                                  newdisease_label=dis)

                    # Normalize by the total number of significant bugs in
                    # that dataset
                    reproducibility = reproducibility/float(coldf.shape[0])
                    results = update_reproducibility_df_lists(
                                  *results,
                                  newvalues=[reproducibility],
                                  newvariables=[dataset],
                                  newmetric_label='rep_dataset_norm',
                                  newdisease_label=dis)
                else:
                    reproducibility = np.nan
                    results = update_reproducibility_df_lists(
                                  *results,
                                  newvalues=[reproducibility],
                                  newvariables=[dataset],
                                  newmetric_label='rep_dataset',
                                  newdisease_label=dis)
                    results = update_reproducibility_df_lists(
                                  *results,
                                  newvalues=[reproducibility],
                                  newvariables=[dataset],
                                  newmetric_label='rep_dataset_norm',
                                  newdisease_label=dis)

            ## Also, if overall is given, calculate the specificity
            # (i.e. how much overlap with the "core" response each dataset has)
            if overall is not None:
                # Note: overall can be a pandas series or one-column dataframe
                overall_healthy = overall[(overall == -1).values].index
                overall_disease = overall[(overall == 1).values].index
                # Note: not looking at overall_mixed bc these *would* be
                # interesting disease-specific bugs! :)

                for dataset in disdf:
                    # Keep just the genera which are significant in the study
                    coldf = disdf[dataset].apply(sigmap)
                    coldf = coldf[coldf != 0]
                    # If there's at least one significant bug in that study
                    if coldf.shape[0] > 0:

                        healthy = coldf[coldf == -1].index
                        disease = coldf[coldf == 1].index

                        healthy_overlap = \
                            len([i for i in healthy if i in overall_healthy])
                        disease_overlap = \
                            len([i for i in disease if i in overall_disease])
                        total_overlap = healthy_overlap + disease_overlap
                        total_sig = len(healthy) + len(disease)
                        total_nonoverlap = total_sig - total_overlap

                        # I don't think I should ever get this error...
                        try:
                            perc_overlap = total_overlap / float(total_sig)
                            perc_nonoverlap = 1.0 - perc_overlap
                        except ZeroDivisionError:
                            perc_overlap = np.nan
                            perc_nonoverlap = np.nan

                    # If nothing was significant, all of these metrics are NaN
                    else:
                         total_overlap = np.nan
                         total_nonoverlap = np.nan
                         total_sig = np.nan
                         perc_overlap = np.nan
                         perc_nonoverlap = np.nan

                    # Update the big results with these new metrics
                    all_overlaps = [total_overlap, total_nonoverlap,
                                    total_sig, perc_overlap,
                                    perc_nonoverlap]
                    all_labels = ['total_overlap', 'total_nonoverlap',
                                  'total_sig', 'perc_overlap',
                                  'perc_nonoverlap']
                    for newvalue, newlabel in zip(all_overlaps, all_labels):
                        results = update_reproducibility_df_lists(
                                      *results,
                                      newvalues=[newvalue],
                                      newvariables=[dataset],
                                      newmetric_label=newlabel,
                                      newdisease_label=dis)

    df_results = pd.DataFrame(data=np.column_stack(results),
                              columns=['value', 'label', 'metric', 'disease'])
    return df_results

def get_dysbiosis_df(dfpvals, qthresh, samplesizes, overall, dfauc):
    """
    Calculates the dysbiosis dataframe from dfpvals, qthresh, samplesizes,
    and overall. Also appends results from classifiers (in dfauc) and returns
    tidy dataframe with *all* the metrics.

    Parameters
    ----------
    dfpvals : pandas dataframe
        genera in rows, datasets in columns, and signed q-values in values.
        Note that 'edd_singh' should be converted to 'cdi_singh' for all
        of these inputs.
    qthresh : float
        significance threshold
    samplesizes : pandas dataframe
        dataframe with datasets in rows, ['total'] sample size column.
    overall : pandas Series
        genera in index and +/-1 in values
        -1 means overall significant in healthy patients
        +1 means overall significant in disease patients
    dfauc : pandas dataframe
        tidy dataframe with 'roc_auc' and 'dataset' columns

    Returns
    -------
    dysbiosis : pandas dataframe
        tidy dataframe with columns 'value', 'label', 'metric', and 'disease'.
    """
    ## Get the tidy dataframe with the "dysbiosis" metrics
    datasets = dfpvals.columns.tolist()
    diseases = list(set([i.split('_')[0] for i in datasets]))
    dysbiosis = get_dysbiosis_metrics(diseases, datasets, dfpvals, qthresh, samplesizes, overall)

    ## Add AUCs to dysbiosis tidy df
    auc_lst = []
    dfauc = dfauc[['dataset', 'roc_auc']].drop_duplicates()
    for d, auc in zip(dfauc['dataset'], dfauc['roc_auc']):
        auc_lst.append([auc, d, 'auc', d.split('_')[0]])
    tidy_auc = pd.DataFrame(data=auc_lst,
                            columns=['value', 'label', 'metric', 'disease'])

    dysbiosis = pd.concat((dysbiosis, tidy_auc))
    dysbiosis['value'] = dysbiosis['value'].astype(float)

    return dysbiosis

p = argparse.ArgumentParser()
p.add_argument('qvalues', help='file with genera in rows, datasets in '
               + 'columns, and signed q-values in values.')
p.add_argument('dataset_info', help='file with datasets in rows, and at '
               + 'least one column "total" with total sample size.')
p.add_argument('overall', help='file with genera in rows, one column named '
               + '"overall", and +/- 1 in values, indicating whether a genus '
               + 'is part of part of the core response.')
p.add_argument('rf_results', help='file with random forest classifier '
               + 'results. Should have columns "dataset" and "roc_auc".')
p.add_argument('dysbiosis_out', help='outfile for tidy file with all of '
               + 'the "dysbiosis metrics".')
p.add_argument('--qthresh', help='significance threshold [default: '
               + ' %(default)s]', default=0.05)
args = p.parse_args()

dfpvals = pd.read_csv(args.qvalues, sep='\t', index_col=0)
samplesizes = pd.read_csv(args.dataset_info, sep='\t', index_col=0)
overall = pd.read_csv(args.overall, sep='\t', index_col=0)
dfauc = pd.read_csv(args.rf_results, sep='\t')
qthresh = args.qthresh

# Need to convert edd_singh to cdi_singh for pattern-matching purposes...
dfpvals = dfpvals.rename(columns={'edd_singh': 'cdi_singh',
    'noncdi_schubert': 'cdi_schubert2'})
samplesizes = samplesizes.rename(index={'edd_singh': 'cdi_singh',
    'noncdi_schubert': 'cdi_schubert2'})
dfauc = dfauc\
    .replace('edd_singh', 'cdi_singh')\
    .replace('noncdi_schubert', 'cdi_schubert2')

dysbiosis = get_dysbiosis_df(dfpvals, qthresh, samplesizes, overall, dfauc)

# Okay switch back to edd_singh
dysbiosis = dysbiosis\
    .replace('cdi_singh', 'edd_singh')\
    .replace('cdi_schubert2', 'noncdi_schubert')

dysbiosis.to_csv(args.dysbiosis_out, sep='\t', index=False)
