#!/usr/bin/env python
"""
This script reads in all of the summary files and writes a few
tables indicating what processing was done to each dataset.
"""

import argparse
import pandas as pd
import numpy as np

# Add this repo to the path
import os, sys
src_dir = os.path.normpath(os.path.join(os.getcwd(), 'src/util'))
sys.path.insert(0, src_dir)
from SummaryParser import SummaryParser
import Formatting as fmt
import FileIO as fio

parser = argparse.ArgumentParser(description='Read and write processing info about the datasets.')
parser.add_argument('yaml_file', help='yaml file with all dataset info')
parser.add_argument('raw_data_dir', help='directory with raw data')
parser.add_argument('proc_info', help='out file with processing info (tex file)')
parser.add_argument('data_info', help='out file with metadata and data sources (tex file)')
args = parser.parse_args()

yamlinfo = fio.read_yaml(args.yaml_file, args.raw_data_dir)

## Get the processing information, as in the yaml file
procinfolst = []
for dataset in yamlinfo:
    sumobj = SummaryParser(yamlinfo[dataset]['summary_file'])
    sumobj.ReadSummaryFile()

    ## What kind of data did we start with?
    try:
        startingfiles = sumobj.attribute_value_16S['RAW_FASTQ_FILE']
        filetype = 'fastq'
        demultiplexed = 'No'
    except:
        try:
            startingfiles = sumobj.attribute_value_16S['RAW_FASTA_FILE']
            filetype = 'fasta'
            demultiplexed = 'No'
        except:
            try:
                startingfiles = sumobj.attribute_value_16S['RAW_FASTQ_FILES']
                filetype = 'fastq'
                demultiplexed = 'Yes'
            except:
                try:
                    startingfiles = sumobj.attribute_value_16S['RAW_FASTA_FILES']
                    filetype = 'fasta'
                    demultiplexed = 'Yes'
                except:
                    raise ValueError('No files specified?')

    if filetype != 'fasta':
        # Barcodes removed?
        barcodes = sumobj.attribute_value_16S['BARCODES_MAP']
        if barcodes == 'None':
            barcodes = 'No'
        else:
            barcodes = 'Yes'

        # Primers trimmed?
        primers = sumobj.attribute_value_16S['PRIMERS_FILE']
        if primers == 'None':
            primers = 'No'
        else:
            primers = 'Yes'

        # Quality filtering or truncating
        try:
            quality = sumobj.attribute_value_16S['QUALITY_TRIM']
            quality_type = '-fastq_truncqual'
        except:
            try:
                quality = sumobj.attribute_value_16S['MAX_ERRORS']
                quality_type = '-fastq_maxee'
            except:
                quality = 25
                quality_type = '-fastq_truncqual'
    else:
        barcodes = np.nan
        primers = np.nan
        quality_type = np.nan
        quality = np.nan

    # Length trim
    length = sumobj.attribute_value_16S['TRIM_LENGTH']

    # Where did the data come from?
    data_source = yamlinfo[dataset]['data_source']
    metadata_source = yamlinfo[dataset]['metadata_source']

    procinfolst.append([dataset, filetype, demultiplexed,
                     barcodes, primers, quality_type, quality, length,
                     data_source, metadata_source])

procinfodf = pd.DataFrame(procinfolst, columns=['dataset',
                                             'data_type', 'demultiplexed',
                                             'barcodes_removed',
                                             'primers_trimmed',
                                             'quality_method',
                                             'quality_cutoff',
                                             'length_trim',
                                             'data_source',
                                             'metadata_source'])

# Sort datasets by dataset ID
procinfodf = procinfodf.sort_values(by='dataset')

# Get human-readable dataset ID
labeldict = fmt.get_labeldict(procinfodf['dataset'])
procinfodf['dataset_label'] = procinfodf['dataset']\
                                  .apply(lambda x: labeldict[x])

procinfodf = procinfodf.replace(np.nan, 'n/a')

## Write the processing info tables
keepcols = ['dataset_label', 'data_type', 'barcodes_removed',
            'primers_trimmed', 'quality_method','quality_cutoff',
            'length_trim']

fmt.write_latex_table(procinfodf[keepcols], args.proc_info)
fmt.write_markdown_table(procinfodf[keepcols],
                         args.proc_info.rsplit('.', 1)[0] + '.md')
procinfodf[keepcols].to_csv(args.proc_info.rsplit('.', 1)[0] + '.txt',
                            sep='\t', index=False)

## Write the data source tables
keepcols = ['dataset_label', 'data_source', 'metadata_source']
fmt.write_latex_table(procinfodf[keepcols], args.data_info)
fmt.write_markdown_table(procinfodf[keepcols],
                         args.data_info.rsplit('.', 1)[0] + '.md')
procinfodf[keepcols].to_csv(args.data_info.rsplit('.', 1)[0] + '.txt',
                            sep='\t', index=False)
