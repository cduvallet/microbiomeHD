#!/usr/bin/env python
"""
This file cleans up the raw OTU and metadata tables and
writes datasetID.otu_table.clean datasetID.metadata.clean.
"""
import argparse
import yaml
import os
import sys
import subprocess

import pandas as pd
import feather
from pyarrow.compat import pdapi

src_dir = os.path.normpath(os.path.join(os.getcwd(), 'src/util'))
sys.path.append(src_dir)
from FileIO import read_yaml

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('raw_data_dir', help='directory with raw results folders')
    p.add_argument('yaml_file', help='yaml file. Should have datasetID '
        + 'as the main key. Can include "condition" to '
        + 'indicate subset of samples to keep.')
    p.add_argument('otu_out', help='clean OTU table output file')

    # Optional args below. Makefile just uses defaults
    p.add_argument('--n-reads-sample', help='minimum reads per sample '
        + ' (default: %(default)s)', default=100)
    p.add_argument('--n-reads-otu', help='minimum reads per OTU '
        + '(default: %(default)s)', default=10)
    p.add_argument('--perc-samples', help='minimum percent of samples an OTU'
        + ' is found in (default: %(default)s)', default=0.01)

    return p.parse_args()

def read_raw_files(otufile, metafile):

    df = pd.read_csv(otufile, sep='\t', index_col=0)
    meta = pd.read_csv(metafile, sep='\t', index_col=0)

    # If either index wasn't read as a string, explicitly do so
    if meta.index.dtype != 'O':
        meta.index = pd.read_csv(y[dataset_id]['metadata_file'], sep='\t', dtype=str).iloc[:,0]

    if df.index.dtype != 'O':
        df.index = pd.read_csv(clean_otu_file, sep='\t', dtype=str).iloc[:,0]

    return df.T, meta


def add_info_to_meta(meta, data, dataset):
    """
    Adds some dataset-specific info to the metadata file.

    Parameters
    ----------
    meta : pandas dataframe
        samples in rows, metadata labels in columns.
    data : dict
        dictionary with specific dataset's entry from the yaml file.
        should have the following keys: 'sequencer', 'region', 'year',
        'disease_label'. These are either in the original yaml file,
        or added by read_yaml()'s default values.
    dataset : str
        dataset ID

    Returns
    -------
    meta : pandas dataframe
        metadata with additional columns
    """

    meta['dataset'] = dataset
    meta['sequencer'] = data['sequencer']
    meta['region'] = data['region']
    meta['year'] = data['year']

    # If the disease label is not DiseaseState, replace it with DiseaseState
    dislabel = data['disease_label']
    if dislabel != 'DiseaseState':
        # If 'DiseaseState' is already in there, delete it
        if 'DiseaseState' in meta.columns:
            print('Removing original DiseaseState column \
                   and replacing it with valuse from {} \
                   column'.format(dislabel))
            meta = meta.drop('DiseaseState', axis=1)
        # Replace dislabel with 'DiseaseState'
        meta.columns = [i if i != dislabel else 'DiseaseState' for i in meta.columns]

    return meta

def clean_up_samples(df, meta, data):
    """
    Cleans up samples in the OTU table and metadata dataframes.
    Keeps only samples which have both metadata and 16S data. If a 'condition'
    was given in the yaml file (as recorded in the data['condition']),
    also keeps only samples which have the specified condition.

    data is a dictionary with {'rawdf': <pandas dataframe>, 'meta': <pandas dataframe>}

    datasetdf is the dataset's dictionary from read_yaml(), i.e. datasets[dataset].
    This function checks if datasetdf has a 'condition' specified.

    returns modified data, which has the same format as input but contains the
    cleaned up dataframes.
    Parameters
    ----------
    df, meta : pandas dataframes
        samples in rows, OTUs/metadata labels in columns.
    data : dict
        dictionary with specific dataset's entry from the yaml file.

    Returns
    -------
    df, meta : pandas dataframes
        dataframes containing only samples which have both
        16S data and metadata. If a subset of metadata was
        specified in 'condition' in the yaml file, only returns
        the subset of samples that meet that condition.
    """

    print('\t\tOriginal: {} samples with 16S, {} samples with metadata.'.format(df.shape[0], meta.shape[0]))

    # If a condition is given, keep only samples with that condition
    try:
        # conditions is a dict with {metadata_column: [keep, conditions]}
        conditions = data['condition']
    except:
        conditions = None
        print('\t\tNo metadata subset specified. Keeping all samples.')

    if conditions is not None:
        for col in conditions:
            print('\t\tKeeping only samples with values {} for metadata ' \
                  'column {}'.format(', '.join([str(i) for i in conditions[col]]), col))
            meta = meta[meta[col].isin(conditions[col])]
            print('\t\t\t{} samples left in metadata'.format(meta.shape[0]))


    # Remove samples which don't have both 16S and metadata
    # if len(meta.index) > len(df.index):
    #     keepsmpls = [i for i in meta.index if i in df.index]
    # else:
    #     keepsmpls = [i for i in df.index if i in meta.index]
    keepsmpls = [i for i in df.index if i in meta.index]

    if len(keepsmpls) != len(df.index) or len(keepsmpls) != len(meta.index):
        print('\t\tDataset has {} samples with 16S, \n\t\t\t{} samples with metadata, \n \
              \t\tbut only {} samples with both'.format(len(df.index), len(meta.index), len(keepsmpls)))
    df = df.loc[keepsmpls]
    meta = meta.loc[keepsmpls]

    return df, meta

def clean_up_tables(df, meta, n_reads_otu, n_reads_sample, perc_samples):
    """
    Cleans up the OTU table and metadata dataframes in data.
    Removes samples with fewer than n_reads_sample reads.
    Removes OTUs with fewer than n_reads_otu reads.
    Removes OTUs which are present in fewer than perc_samples*100 percent of samples.
    Removes empty samples/OTUs.
    """

    # Remove samples with fewer than n_reads reads.
    df = remove_shallow_smpls(df, n_reads_sample)

    # Remove OTUs with fewer than 10 reads
    old = df.shape[1]
    df = remove_shallow_otus(df, n_reads=n_reads_otu)
    new = df.shape[1]
    if new < old:
        print('\t\tOf {} original OTUs, {} have more than {} reads'.format(old, new, n_reads_otu))

    # Remove OTUs which are present in fewer than perc_samples of samples.
    old = df.shape[1]
    df = remove_shallow_otus(df, perc_samples=perc_samples)
    new = df.shape[1]
    if new < old:
        print('\t\tOf {} original OTUs, {} are present \n\t\t\tin more than {}% of samples'.format(old, new, perc_samples*100))

    # Remove any samples which now have fewer than n_reads
    df = remove_shallow_smpls(df, n_reads_sample)

    # Double check that both metadata and OTU table have same samples
    # (after we've filtered out samples based on reads, etc)
    # if len(meta.index) > len(df.index):
    #     keepsmpls = [i for i in meta.index if i in df.index]
    # else:
    #     keepsmpls = [i for i in df.index if i in meta.index]
    keepsmpls = [i for i in df.index if i in meta.index]

    if len(keepsmpls) != len(df.index) or len(keepsmpls) != len(meta.index):
        print('\t\tAfter some cleaning, dataset has {} samples with 16S, \n \
              \t\t{} samples with metadata, but only {} samples with both'.format(
              len(df.index), len(meta.index), len(keepsmpls)))
        df = df.loc[keepsmpls]
        meta = meta.loc[keepsmpls]

    # Remove empty metadata columns (just in case...)
    meta = meta.dropna(how='all', axis=1)

    return df, meta

def remove_shallow_smpls(df, n_reads):
    """
    Removes samples with fewer than n_reads from dataframe df.

    Parameters
    -----------
    df : pandas dataframe
        samples are in rows, OTUs are in columns
    n_reads : int
        minimum number of reads per sample for sample to be kept
    """

    total_reads = df.sum(axis=1)
    shallow_smpls = [smpl for smpl in total_reads.index \
                     if total_reads.loc[smpl] <= n_reads]
    df = df.drop(shallow_smpls)

    return df

def remove_shallow_otus(df, perc_samples=None, n_reads=None):
    """
    Removes OTUs which are present in fewer than 100*perc_samples percent of
    samples OR which have fewer than n_reads reads.

    Parameters
    ----------
    df : pandas dataframe
        Samples are in rows. OTUs are in columns.
    perc_samples : float
        min percent of samples that an OTU must be present in to not
        be thrown out.
    n_reads : int
        min number of reads an OTU must have in df to not be thrown
        out.

    Either perc_samples or n_reads must be specified. If both are specified,
    the perc_samples filtering is done first and then OTUs with fewer than
    n_reads total are thrown out.

    """
    if perc_samples is not None:
        presencemap = lambda x: 1 if x else 0
        otus_perc_present = df.applymap(presencemap).sum() / df.shape[0]
        keepotus = list(
            otus_perc_present[otus_perc_present > perc_samples].index)
        df = df[keepotus]

    if n_reads is not None:
        # Removes any OTUs with fewer than n_reads from the raw and abun dfs
        # samples are in rows and OTUs are in columns
        total_reads = df.sum(axis=0)
        shallow_col_indices = [i for i in range(len(total_reads.index)) \
                               if total_reads.iloc[i] < n_reads]
        shallow_otus = df.columns[shallow_col_indices]
        df = df.drop(shallow_otus, axis=1)

    return df

def fix_ob_zhu(meta):
    """
    Fix the DiseaseState labels for case patients in the ob_zhu data.
    """
    meta['DiseaseState'] = meta['DiseaseState'].replace('nonNASH-OB', 'OB')
    meta['DiseaseState'] = meta['DiseaseState'].replace('NASH', 'OB-NASH')
    return meta

def fix_cdi_schubert(meta):
    """
    Fix the DiseaseState labels for case patients in cdi_schubert
    """
    meta['DiseaseState'] = meta['DiseaseState']\
        .replace('nonCDI', 'ignore-nonCDI')
    return meta

def fix_noncdi_schubert(meta):
    meta['DiseaseState'] = meta['DiseaseState'].replace('CDI', 'ignore-CDI')
    return meta

if __name__ == "__main__":

    args = parse_args()

    dataset_id = args.otu_out.split('/')[-1].split('.')[0]
    y = read_yaml(args.yaml_file, args.raw_data_dir)

    df, meta = read_raw_files(y[dataset_id]['otu_table'],
                              y[dataset_id]['metadata_file'])

    ## Add some study-wise metadata, like sequencer and region
    meta = add_info_to_meta(meta, y[dataset_id], dataset_id)

    df, meta = clean_up_samples(df, meta, y[dataset_id])
    df, meta = clean_up_tables(df, meta, args.n_reads_otu, args.n_reads_sample, args.perc_samples)

    ## Add sequencing depth to metadata
    meta['total_reads'] = df.sum(axis=1)

    if dataset_id == 'ob_zhu':
        # Turn the nonNASH-OB into OB and the NASH into OB-NASH,
        # so that they're not recognized as cases in downstream analyses
        meta = fix_ob_zhu(meta)
    elif dataset_id == 'cdi_schubert':
        meta = fix_cdi_schubert(meta)
    elif dataset_id == 'noncdi_schubert':
        meta = fix_noncdi_schubert(meta)

    # Reset indices to write as feather format
    df = df.reset_index()
    feather.write_dataframe(df, args.otu_out)

    meta_out = args.otu_out.split('.otu_table.clean.feather')[0] + '.metadata.clean.feather'
    meta = meta.reset_index()

    # Feather doesn't support writing Object column types with 'mixed' inferred
    # dtype OR non-unicode or non-string inferred dtype.
    # Need to convert any Object columns to their inferred dtype.
    # https://github.com/wesm/feather/blob/master/python/feather/api.py#L42
    for i, name in enumerate(meta.columns):
        col = meta.iloc[:, i]

        if pdapi.is_object_dtype(col):
            inferred_type = pd.lib.infer_dtype(col)
            if inferred_type == "boolean":
                meta.iloc[:, i] = meta.iloc[:, i].astype(bool)

    feather.write_dataframe(meta, meta_out)
