#!/usr/bin/env python

"""
This script reads in all the OTU tables and metadata and the
yaml file and writes a file with some information about each
dataset.
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
from SummaryParser import SummaryParser

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

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('yaml_file', help='yaml file with all dataset info')
    parser.add_argument('raw_data_dir', help='directory with raw data')
    parser.add_argument('clean_data_dir', help='directory with clean OTU tables and metadata')
    parser.add_argument('dataset_info', help='out file with dataset info')
    parser.add_argument('proc_info', help='out file with processing info')
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
    # Save
    stats.to_csv(args.dataset_info, sep='\t', index=False)

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

    procinfodf.to_csv(args.proc_info, sep='\t', index=False)
