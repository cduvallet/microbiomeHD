#!/usr/bin/env
"""
This scripts writes a markdown table describing all (?) the datasets in
Zenodo, for the Zenodo description.
"""
import requests
import argparse
import yaml
import pandas as pd

# Add util repo to the path
import os, sys
src_dir = os.path.normpath(os.path.join(os.getcwd(), 'src/util'))
sys.path.insert(0, src_dir)
import Formatting as fmt


p = argparse.ArgumentParser()
p.add_argument('token_file', help='Zenodo access token file')
p.add_argument('yaml_file', help='path to yaml file, which has tar file '
    + '--> datasetID')
p.add_argument('out', help='out file')
args = p.parse_args()

# Get the
with open(args.token_file, 'r') as f:
    token = f.readline().strip()

r = requests.get('https://zenodo.org/api/deposit/depositions/1146764/files',
                 params={'access_token': token})

names = [i['filename'].split('.tar.gz')[0]
         for i in r.json() if i['filename'].endswith('.tar.gz')]

# Read in yaml file
with open(args.yaml_file, 'r') as f:
    y = yaml.load(f)

folder2id = {y[i]['folder']: i for i in y}

# folder, dataset ID, controls, n controls, disease, n cases, paper
table_list = []
for f in names:
    fname = f + '.tar.gz'
    datasetid = folder2id[f]
    paper = y[datasetid]['paper']
    # This returns a dict, which I need to convert to a str with <br>'s
    samples = y[datasetid]['sample_size']
    samples = ', '.join(
        ["{}: {}".format(i, str(samples[i])) for i in samples])
    table_list.append([fname, datasetid, samples, paper])

with open(args.out, 'w') as f:
    for line in table_list:
        f.write('* **{}** (*{}*): {}\n\t* {}\n'.format(*line))
