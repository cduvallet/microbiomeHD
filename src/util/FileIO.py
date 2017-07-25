#!/usr/bin/env python
"""
Functions to interface with the OTU and metadata files.
"""

import os, sys
import yaml
import pandas as pd
import feather

# Add this repo to the path
src_dir = os.path.normpath(os.path.join(os.getcwd(), 'src/util'))
sys.path.insert(0, src_dir)
from util import raw2abun

def read_yaml(yamlfile, batch_data_dir):
    """
    Reads in a yaml file with {dataset_id: {
                                folder: <folder name>
                                region: <16S Region>
                                sequencer: <DNA sequencer>
                                condition: <condition_dict>},
                                disease_label: <label>,
                                table_type: <'classic' or 'normal'>,
                              dataset2: {...}, ...}

    Returns a dict with {dataset_id: {
                                otu_table: <path to otu table file>,
                                meta_file: <path to metadata file>,
                                summary_file: <path to summary_file.txt>,
                                sequencer: <DNA sequencer>,
                                region: <16S region>,
                                disease_label: 'label',
                                table_type: <'classic' or 'normal'>,
                                condition: <condition_dict, if given>, ...}
    yaml file can have 'otu_table' and 'metadata_file' keys indicating full paths
        to the respective files.

        Otherwise, it can have a 'folder' key which indicates the results_folder
        name in batch_data_dir.
        If a 'folder' is given, otu_table and metadata files are assumed to be
                              folder/RDP/<folder minus '_results'>.otu_table.100.denovo.rdp_assigned
                              folder/<folder minus '_results'>.metadata.txt

        'table_type' key indicates whether samples or OTUs are in rows.
            'normal' means that OTUs are in columns and samples are in rows
            'classic' means that OTUs are in rows and samples are in columns
        If not given, defaults to 'classic'
    """

    with open(yamlfile, 'r') as f:
        datasets = yaml.load(f)

    for dataset in datasets:
        # Grab the sub-dict with just that dataset's info
        data = datasets[dataset]

        # Get OTU table file, if it's not already specified
        if 'otu_table' not in data:
            try:
                folder = data['folder']
                folderpath = os.path.join(batch_data_dir, folder)
                # folderpath/RDP/datasetID.otu_table.100.denovo.rdp_assigned
                datasets[dataset]['otu_table'] = \
                    os.path.realpath(
                        os.path.join(os.path.join(folderpath, 'RDP'),
                                     folder.rsplit('_', 1)[0]
                                     + '.otu_table.100.denovo.rdp_assigned'))
            except KeyError:
                raise ValueError('No OTU table or results folder specified for dataset {}'.format(dataset))

        # Get metadata file, if it's not already specified
        if 'metadata_file' not in data:
            try:
                folder = data['folder']
                folderpath = os.path.join(batch_data_dir, folder)
                # folderpath/datasetID.metadata.txt
                datasets[dataset]['metadata_file'] = \
                    os.path.realpath(
                        os.path.join(folderpath,
                                     folder.rsplit('_',1)[0] + '.metadata.txt'))
            except KeyError:
                raise ValueError('No metadata file or results folder specified for dataset {}'.format(dataset))

        # Get summary_file.txt
        if 'summary_file' not in data:
            try:
                folder = data['folder']
                folderpath = os.path.join(batch_data_dir, folder)
                # folderpath/summary_file.txt
                datasets[dataset]['summary_file'] = \
                        os.path.realpath(
                            os.path.join(folderpath, 'summary_file.txt'))
            except KeyError:
                print('No summary file or results folder specified for dataset {}'.format(dataset))

        # If 16S region is not specified, add it to dict as 'unk'
        if 'region' not in data:
            datasets[dataset]['region'] = 'unk'

        # If sequencer is not specified, add it to dict as 'unk'
        if 'sequencer' not in data:
            datasets[dataset]['sequencer'] = 'unk'

        # diseaselabel is not specified, assumes that it's DiseaseState
        if 'disease_label' not in data:
            datasets[dataset]['disease_label'] = 'DiseaseState'

        if 'table_type' not in data:
            datasets[dataset]['table_type'] = 'classic'

        if 'year' not in data:
            datasets[dataset]['year'] = 'unk'

    return datasets

def read_dataset_files(datasetid, clean_folder):
    """
    Reads the OTU table and metadata files for datasetid in clean_folder.
    """
    fnotu = datasetid + '.otu_table.clean.feather'
    fnmeta = datasetid + '.metadata.clean.feather'

    df = feather.read_dataframe(os.path.join(clean_folder, fnotu))
    # Feather format does not support index names, first column has index
    df.index = df.iloc[:,0]
    df = df.iloc[:, 1:]

    meta = feather.read_dataframe(os.path.join(clean_folder, fnmeta))
    meta.index = meta.iloc[:, 0]
    meta = meta.iloc[:, 1:]

    ## Make sure sample names are strings
    if df.index.dtype != 'O':
        df.index = pd.read_csv(os.path.join(clean_folder, fnotu), sep='\t', dtype=str).iloc[:,0]

    if meta.index.dtype != 'O':
        meta.index = pd.read_csv(os.path.join(clean_folder, fnmeta), sep='\t', dtype=str).iloc[:,0]

    return df, meta

def get_classes(meta, dataset=''):
    """
    Returns classes_list for supervised comparison. List of accepted controls
    and diseases is hard-coded in this function.

    Looks in meta['DiseaseState'] to return
    [[control label(s)], [disease label(s)]]

    Parameters
    ----------
    meta : pandas DataFrame
        metadata with 'DiseaseState' column
    dataset : str
        dataset ID, used in ValueError message, if one of the classes has
        no labels
    noncdi : bool
        whether to include nonCDI in the diseases list

    Returns
    -------
    classes_list : list
        [list of control lables, list of case labels]
    """
    controls = ['H', 'nonIBD']
    diseases = ['ASD', 'CD', 'CDI', 'nonCDI', 'CIRR', 'CRC', 'EDD', 'HIV',
                'MHE', 'NASH', 'OB', 'PAR', 'PSA', 'RA', 'T1D', 'T2D',
                 'UC']

    labels = list(set(meta['DiseaseState']))

    classes_list = [[i for i in labels if i in controls],
                    [j for j in labels if j in diseases]]

    if len(classes_list[0]) == 0 or len(classes_list[1]) == 0:
        raise ValueError('Something wrong with ' + dataset + ' DiseaseState.')

    return classes_list

# def get_separate_case_classes(meta, dataset=''):
#     """
#     Returns list of case labels for diseases which need expanded analyses
#     that consider the heterogeneity of case patients.
#
#     This is mostly a stand-in for a place to hard-code in some diseases,
#     which are different than the above get_classes() function (because, for
#     example, we don't want to include nonCDI cases in that function but we
#     do here).
#
#     Parameters
#     ----------
#     meta : pandas DataFrame
#         metadata with 'DiseaseState' column
#     dataset : str
#         dataset ID, used in ValueError message, if one of the classes has
#         no labels
#
#     Returns
#     -------
#     classes_list : list
#         [list of control lables, list of case labels]
#     """
#     controls = ['H', 'nonIBD']
#     diseases = ['nonCDI', 'UC', 'CD', 'CIRR', 'MHE', 'RA', 'PSA', 'EDD', 'CDI']
#
#     labels = list(set(meta['DiseaseState']))
#
#     classes_list = [[i for i in labels if i in controls],
#                     [j for j in labels if j in diseases]]
#
#     if len(classes_list[0]) == 0 or len(classes_list[1]) == 0:
#         raise ValueError('Something wrong with ' + dataset + ' DiseaseState.')
#
#     return classes_list

def get_samples(meta, classes_list):
    """
    Returns [[healthy samples], [disease samples]] from metadata corresponding to
    lables in classes_list

    Parameters
    ----------
    meta          pandas df with 'DiseaseState' column
    classes_list  list, [['DiseaseState' label for control samples], ['DiseaseState' label for case samples]]

    Returns
    -------
    list, [[control patients], [diseases patients]]
    """

    h = list(meta[meta['DiseaseState'].isin(classes_list[0])].index)
    dis = list(meta[meta['DiseaseState'].isin(classes_list[1])].index)

    return [h, dis]

def get_dataset_ids(clean_folder):
    """
    Gets the list of dataset_ids present in clean_folder.
    """
    files = os.listdir(clean_folder)

    datasets = list(set([i.split('.')[0] for i in files]))

    return [d for d in datasets if d + '.otu_table.clean.feather' in files and d + '.metadata.clean.feather' in files]

def read_dfdict_data(datadir, subset=None):
    """
    Read in all df's, metadata, dis_smpls, H_smpls, and classes_list for all
    datasets in datadir.

    Parameters
    ----------
    datadir : str
        path to directory with all of the dataset.otu_table.clean
        and dataset.metadata.clean files
    subset : str
        path to file with subset of datasets to read. Each dataset ID should
        be on one line, and should have a corresonding fnotu and fnmeta,
        as defined in read_dataset_files() function.

    Returns
    -------
    dfdict : dict
        {dataset: {'df': df, 'meta': meta, 'dis_smpls': list_of_dis_smpls,
                   'H_smpls': list_of_H_smpls, 'classes': classes_list}}
    """
    print('Reading datasets...')
    # Initialize dict to store all dataframes
    dfdict = {}

    # If subset of datasets are given, read only those
    if subset is not None:
        with open(subset, 'r') as f:
            datasetids = f.read().splitlines()
    else:
        datasetids = get_dataset_ids(datadir)

    # Read each dataset and convert to relative abundance
    for dataset in datasetids:
        print(dataset),
        ## Read dataset
        df, meta = read_dataset_files(dataset, datadir)
        df = raw2abun(df)

        ## Get case and control samples
        classes_list = get_classes(meta)
        if len(classes_list[0]) == 0 or len(classes_list[1]) == 0:
            raise ValueError('Something wrong with ' + dataset + ' metadata.')
        H_smpls, dis_smpls = get_samples(meta, classes_list)

        dfdict.update({dataset: {'df': df, 'meta': meta, 'dis_smpls': dis_smpls, 'H_smpls': H_smpls, 'classes': classes_list}})
    print('\nReading datasets... Finished.')
    return dfdict
