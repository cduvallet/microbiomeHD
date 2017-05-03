#!/usr/bin/env python
"""
This script makes Table 1 and Supplementary Table 1, which have the basic
dataset info like number of samples, citation, and basic sequencing info.
"""

import os
import sys
import pandas as pd
import numpy as np
import argparse

# Add this repo to the path
src_dir = os.path.normpath(os.path.join(os.getcwd(), 'src/util'))
sys.path.append(src_dir)
import FileIO as fio
import Formatting as fmt

def get_citation(dataset):
    """
    Get citation for each datasets, as in my refs.bib file
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
         'ra_scher': 'ra-littman',
         't1d_mejialeon': 't1d-mejia'}
    if dataset in d:
        return d[dataset]
    else:
        return '-'.join(dataset.split('_'))

parser = argparse.ArgumentParser(description='Read and write basic info about the datasets.')
parser.add_argument('yaml_file', help='yaml file with all dataset info')
parser.add_argument('raw_data_dir', help='directory with raw data')
parser.add_argument('clean_data_dir', help='directory with clean OTU tables and metadata')
parser.add_argument('dataset_info', help='out file with the tab-delimited, '
                    + 'unformatted information.')
parser.add_argument('main_table', help='out file with table for main text (tex file)')
parser.add_argument('supp_table', help='out file with supplementary table (tex file)')
args = parser.parse_args()

yamlinfo = fio.read_yaml(args.yaml_file, args.raw_data_dir)

## Get the dataset information, like sample size and read depth
statslst = []
for dataset in yamlinfo:
    print(dataset),
    df, meta = fio.read_dataset_files(dataset, args.clean_data_dir)

    classes_list = fio.get_classes(meta)
    [h_smpls, dis_smpls] = fio.get_samples(meta, classes_list)

    # String with the types of patients in control/case categories
    controls = ', '.join(classes_list[0])
    cases = ', '.join(classes_list[1])

    sequencer = yamlinfo[dataset]['sequencer']
    year = yamlinfo[dataset]['year']
    region = yamlinfo[dataset]['region']

    # Citation for Latex string
    citation = "\cite{" + get_citation(dataset) + "}"

    statslst.append([dataset,
                     len(h_smpls), len(dis_smpls),
                     len(h_smpls) + len(dis_smpls),
                     controls, cases,
                     df.sum(axis=1).min(),
                     df.sum(axis=1).max(),
                     df.sum(axis=1).median(),
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

## Format and save Table 1 latex and markdown files

# Re-order dataframe to match order in Figure 1
stats['dataset'] = stats['dataset'].apply(lambda x: x if x != 'edd_singh' else 'cdi_singh')
_, dataset_order = fmt.get_dataset_order(stats)
stats['dataset'] = stats['dataset'].astype('category')
stats['dataset'].cat.set_categories(dataset_order, inplace=True)
stats = stats.sort_values(by='dataset')

labeldict = fmt.get_labeldict(stats['dataset'])
stats['dataset_label'] = stats['dataset'].apply(lambda x: labeldict[x])

table1 = stats[['dataset_label', 'N_ctrl', 'controls', 'N_dis',
                'cases', 'citation']]
fmt.write_latex_table(table1, args.main_table)
fmt.write_markdown_table(table1, args.main_table.rsplit('.', 1)[0] + '.md')
table1.to_csv(args.main_table.rsplit('.', 1)[0] + '.txt', sep='\t', index=False)

## Format and save Table 2 latex and markdown files
# Re-order alphabetically by dataset ID (i.e. disease first,
# then alpha by author)
stats['dataset'] = stats['dataset'].astype(str)
stats = stats.sort_values(by='dataset')
table2 = stats[['dataset_label', 'year', 'N_ctrl', 'controls',
                'N_dis', 'cases', 'med_reads', 'sequencer',
                'region', 'citation']]
fmt.write_latex_table(table2, args.supp_table)
fmt.write_markdown_table(table2, args.supp_table.rsplit('.',1)[0] + '.md')
table2.to_csv(args.supp_table.rsplit('.',1)[0] + '.txt', sep='\t', index=False)
