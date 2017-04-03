#!/usr/bin/env python

"""
This script cleans up the raw OTU tables from processing.
"""

import argparse
import yaml
import os
import sys
import subprocess

import pandas as pd

src_dir = os.path.normpath(os.path.join(os.getcwd(), 'src/util'))
sys.path.append(src_dir)
from FileIO import read_yaml
#
# def get_tar_file(data):
#     """
#     Return the raw_data tar file or directory with the .
#
#     Parameters
#     ----------
#     data : dict
#         dictionary corresponding to the yaml info for the dataset
#         of interest. Should either have a 'tar_file' key which
#         directly gives the tar file, or a 'folder' key which
#         is assumed to be the prefix for the tar file
#
#     Returns
#     -------
#     tar_file : str
#         the filename with the "raw" OTU tables for this dataset
#     """
#
#     if 'tar_file' in data:
#         return data['tar_file']
#     elif 'folder' in data:
#         return data['folder'] + '.tar.gz'
#     else:
#         raise ValueError("No raw data specified. Make sure your yaml"
#                          + " file has a 'tar_file' or 'folder' key.")

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

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('raw_data_dir', help='directory with raw OTU tables')
    parser.add_argument('yaml_file', help='yaml file. Should have datasetID '
                        + 'as the main key, and at least one of the sub-keys '
                        + '"tar file" (that gives the name of the .tar.gz file '
                        + 'with the corresponding raw OTU tables) or '
                        + '"folder" (that is converted to folder.tar.gz '
                        + 'and is in raw_data_dir.)')
    parser.add_argument('clean_file', help='clean OTU table name. Should have'
                        + ' format: datasetID.otu_table.clean')
    # Optional args below. Makefile just uses defaults
    parser.add_argument('--sample-min', help='minimum reads per sample '
                        + ' (default: %(default)s)', default=100)
    parser.add_argument('--otu-min', help='minimum reads per OTU '
                        + '(default: %(default)s)', default=10)
    parser.add_argument('--otu-perc', help='minimum percent of samples an OTU'
                        + ' is found in (default: %(default)s)', default=0.01)

    args = parser.parse_args()

    dataset_id = args.clean_file.split('/')[-1].split('.')[0]
    y = read_yaml(args.yaml_file, args.raw_data_dir)

    # # Get the tar file and unzip it
    # tar_file = get_tar_file(y[dataset_id])
    # tar_file = os.path.join(args.raw_data_dir, tar_file)
    # cmd = 'tar -C {} -zxvkf {}'.format(args.raw_data_dir, tar_file)
    # subprocess.call(cmd, shell=True)

    # Read the OTU table
    df = pd.read_csv(y[dataset_id]['otu_table'], sep='\t', index_col=0).T

    ## Remove samples and OTUs based on criteria
    # Remove samples with fewer than 100 reads.
    df = remove_shallow_smpls(df, n_reads=args.sample_min)

    # Remove OTUs in fewer than 1% of samples and with fewer than 10 reads
    df = remove_shallow_otus(df, perc_samples=args.otu_perc,
                             n_reads=args.otu_min)

    # Remove any samples which now have fewer than n_reads
    df = remove_shallow_smpls(df, n_reads=args.sample_min)

    # Write the OTU table
    df.to_csv(args.clean_file, sep='\t')

    print(dataset_id, y[dataset_id]['otu_table'], args.clean_file)
