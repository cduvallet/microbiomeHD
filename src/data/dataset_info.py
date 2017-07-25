#!/usr/bin/env python
"""
This script reads in the yaml file and datasets and outputs basic
information about each dataset.
"""

import os
import sys
import pandas as pd
import argparse

# Add this repo to the path
src_dir = os.path.normpath(os.path.join(os.getcwd(), 'src/util'))
sys.path.insert(0, src_dir)
import FileIO as fio

def get_citation(dataset):
    """
    Get citation for each datasets, as in my refs.bib file.
    """
    d = {'asd_kang': 'asd-kb',
         'crc_chen': 'crc-xiang',
         'crc_wang': 'crc-zhao',
         'hiv_lozupone': 'lozupone2013alterations',
         'hiv_noguerajulian': 'noguera2016gut',
         'ibd_morgan': 'ibd-hut',
         'ibd_willing': 'ibd-engstrand',
         'nash_wong': 'nash-chan',
         'nash_zhu': 'nash-baker',
         'ob_turnbaugh': 'ob-gordon',
         'ob_zhu': 'nash-baker',
         'par_scheperjans': 'par-schep',
         'art_scher': 'ra-littman',
         't1d_mejialeon': 't1d-mejia',
         'cd_gevers': 'ibd-gevers',
         'cd_morgan': 'ibd-hut',
         'uc_morgan': 'ibd-hut',
         'cd_papa': 'ibd-papa',
         'uc_papa': 'ibd-papa',
         'cd_willing': 'ibd-engstrand',
         'uc_willing': 'ibd-engstrand',
         'liv_zhang': 'mhe-zhang',
         'cirr_zhang': 'mhe-zhang',
         'ra_scher': 'ra-littman',
         'psa_scher': 'ra-littman',
         'noncdi_schubert': 'cdi-schubert',
         'cdi_schubert2': 'cdi-schubert'}
    if dataset in d:
        return d[dataset]
    else:
        return '-'.join(dataset.split('_'))


parser = argparse.ArgumentParser(description='Read and write basic info '
    + 'about the datasets.')
parser.add_argument('yaml_file', help='yaml file with all dataset info')
parser.add_argument('raw_data_dir', help='directory with raw data')
parser.add_argument('clean_data_dir', help='directory with clean OTU tables '
    + 'and metadata')
parser.add_argument('dataset_info', help='out file with the tab-delimited, '
    + 'unformatted information.')
parser.add_argument('--subset', help='file with list of dataset IDs that '
    + 'should be analyzed, if you want to analyze only a subset of the '
    + 'datasets in clean_data_dir. Dataset IDs should match IDs in file '
    + 'names, and be one per line.', default=None)
parser.add_argument('--split-cases', help='flag to provide sample size info '
    + 'for each case group separately', action='store_true')
args = parser.parse_args()

yamlinfo = fio.read_yaml(args.yaml_file, args.raw_data_dir)

if args.subset is not None:
    with open(args.subset, 'r') as f:
        datasetids = f.read().splitlines()
else:
    datasetids = yamlinfo.keys()

## Get the dataset information, like sample size and read depth
statslst = []
for dataset in datasetids:
    print(dataset),

    # Dataset metadata
    sequencer = yamlinfo[dataset]['sequencer']
    year = yamlinfo[dataset]['year']
    region = yamlinfo[dataset]['region']
    # Citation for Latex string
    citation = "\cite{" + get_citation(dataset) + "}"

    # Data-dependent info
    df, meta = fio.read_dataset_files(dataset, args.clean_data_dir)
    if args.split_cases:
        classes_list = fio.get_classes(meta)

        for dis_label in classes_list[1]:
            # old dataset = ibd_alm, new dataset = uc_alm
            newdataset = dis_label.lower() + '_' + dataset.split('_')[1]

            # Get samples
            sub_list = [classes_list[0], [dis_label]]
            h_smpls, dis_smpls = fio.get_samples(meta, sub_list)
            all_smpls = h_smpls + dis_smpls

            # String with the types of patients in control/case categories
            controls = ', '.join(sub_list[0])
            cases = ', '.join(sub_list[1])

            statslst.append([newdataset,
                             len(h_smpls), len(dis_smpls),
                             len(all_smpls),
                             controls, cases,
                             df.loc[all_smpls].sum(axis=1).min(),
                             df.loc[all_smpls].sum(axis=1).max(),
                             df.loc[all_smpls].sum(axis=1).median(),
                             sequencer, region,
                             year, citation])

    else:
        classes_list = fio.get_classes(meta)
        h_smpls, dis_smpls = fio.get_samples(meta, classes_list)
        all_smpls = h_smpls + dis_smpls

        # String with the types of patients in control/case categories
        controls = ', '.join(classes_list[0])
        cases = ', '.join(classes_list[1])

        statslst.append([dataset,
                         len(h_smpls), len(dis_smpls),
                         len(all_smpls),
                         controls, cases,
                         df.loc[all_smpls].sum(axis=1).min(),
                         df.loc[all_smpls].sum(axis=1).max(),
                         df.loc[all_smpls].sum(axis=1).median(),
                         sequencer, region,
                         year, citation])

# Convert to dataframe
stats = pd.DataFrame(data=statslst,
                     columns=['dataset', 'N_ctrl', 'N_dis',
                              'total', 'controls', 'cases',
                              'min_reads', 'max_reads',
                              'med_reads', 'sequencer',
                              'region', 'year', 'citation'])
# Save raw data
stats.to_csv(args.dataset_info, sep='\t', index=False)
