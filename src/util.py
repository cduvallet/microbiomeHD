#!/usr/bin/env python

"""
Useful functions to be used through data processing and analysis code.
"""

import os
import yaml

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
