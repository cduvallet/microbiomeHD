#!/usr/bin/env python

"""
This script cleans up the raw OTU tables from processing.
"""

import argparse
import yaml
import os
import sys

import pandas as pd

src_dir = os.path.normpath(os.path.join(os.getcwd(), 'src/util'))
sys.path.append(src_dir)
from FileIO import read_yaml

def clean_up_samples_from_meta(df, meta, data):
    """
    Cleans up samples in the OTU table and metadata dataframes.
    Keeps only samples which have both metadata and 16S data. If a 'condition'
    was given in the yaml file (as recorded in the data dictonary),
    also keeps only samples which have the specified condition.

    Parameters
    ----------
    df : pandas dataframe
        OTU table with samples in rows, OTUs in columns
    meta : pandas dataframe
        OTU table with samples in rows, metadata labels in columns
    data : dict
        dictionary with specific dataset's entry in the yaml file,
        i.e. yaml_dict[dataset]

    Returns
    -------
    df, meta : pandas dataframes
        modified df and meta with the correct samples
    """

    print('\t\tOriginal: {} samples with 16S, {} samples with metadata.'.format(df.shape[0], meta.shape[0]))

    # If a condition is given, keep only samples with that condition
    try:
        conditions = data['condition']    # conditions is a dict with {metadata_column: [keep, conditions]}
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
    if len(meta.index) > len(df.index):
        keepsmpls = [i for i in meta.index if i in df.index]
    else:
        keepsmpls = [i for i in df.index if i in meta.index]

    if len(keepsmpls) != len(df.index) or len(keepsmpls) != len(meta.index):
        print('\t\tDataset has {} samples with 16S, \n\t\t\t{} samples with metadata, \n \
              \t\tbut only {} samples with both'.format(len(df.index), len(meta.index), len(keepsmpls)))
    df = df.loc[keepsmpls]
    meta = meta.loc[keepsmpls]

    return df, meta

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

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('raw_data_dir', help='directory with raw metadata')
    parser.add_argument('yaml_file', help='yaml file. Should have datasetID '
                        + 'as the main key. Can include "condition" to '
                        + 'indicate subset of samples to keep.')
    parser.add_argument('clean_file', help='clean metadata file name. '
                        + 'Should have format: datasetID.metadata.clean')

    args = parser.parse_args()

    dataset_id = args.clean_file.split('/')[-1].split('.')[0]
    y = read_yaml(args.yaml_file, args.raw_data_dir)

    # Read the raw metadata and OTU table files
    meta = pd.read_csv(y[dataset_id]['metadata_file'], sep='\t',
                       index_col=0)
    df = pd.read_csv(y[dataset_id]['otu_table'], sep='\t', index_col=0).T

    # Keep samples with both 16S and metadata, and which follow any
    # conditions specified in the yaml file
    df, meta = clean_up_samples_from_meta(df, meta, y[dataset_id])

    # Add manually-defined parameters to the metadata file
    meta = add_info_to_meta(meta, y[dataset_id], dataset_id)

    # Add sequencing depth to metadata
    meta['total_reads'] = df.sum(axis=1)

    # Write both metadata and OTU table
    meta.to_csv(args.clean_file, sep='\t')
    otu_file = args.clean_file.split('.metadata.clean')[0] + '.otu_table.clean'
    df.to_csv(otu_file, sep='\t')
