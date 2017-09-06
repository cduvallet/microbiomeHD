#!/usr/bin/env python
"""
This script calculates the expected number of core genera under the null
model.
"""
import argparse
import pandas as pd
import numpy as np

# Add this repo to the path
import os, sys
src_dir = os.path.normpath(os.path.join(os.getcwd(), 'src/analysis'))
sys.path.insert(0, src_dir)
src_dir = os.path.normpath(os.path.join(os.getcwd(), 'src/util'))
sys.path.insert(0, src_dir)
from meta_analyze import count_sig, cross_disease_meta_analysis
from util import shuffle_col

parser = argparse.ArgumentParser()
parser.add_argument('qvalues', help='file with qvalues; genera in rows, '
    + ' datasets in columns')
parser.add_argument('qthresh', help='significance threshold', default=0.05,
    type=float)
parser.add_argument('out', help='file to write tidy results to')
parser.add_argument('--n_diseases', help='number of diseases to use in '
    + ' calculating "core" genera', default=2, type=int)
parser.add_argument('--exclude-nonhealthy', help='flag to exclude '
    + 'studies without healthy controls and hiv_lozupone from the '
    + 'overall cross-disease meta-analysis', action='store_true')
parser.add_argument('--reps', help='number of repetitions to build null '
    + '[default: %(default)s]', default=1000, type=int)

args = parser.parse_args()

qvals = pd.read_csv(args.qvalues, sep='\t', index_col=0)

if args.exclude_nonhealthy:
    to_exclude = ['ibd_papa', 'ibd_gevers', 'hiv_lozupone']
    qvals = qvals.drop(to_exclude, axis=1)
    meta_counts = count_sig(qvals, args.qthresh)

# Placeholder for removing datasets to exclude... should actually probably
# read these in and include their removal in count_sig call...

# For each repetition, shuffle labels, count number of sig bugs
results = []
for i in range(args.reps):
    print(i),
    newq = qvals.copy().apply(shuffle_col)
    counts = count_sig(newq, args.qthresh)
    overall_df = cross_disease_meta_analysis(counts, args.n_diseases)
    for c, n in zip(['health', 'mixed', 'disease'], [-1, 0, 1]):
        try:
            results.append([i, c, overall_df.groupby('overall').size()[n]])
        except KeyError:
            results.append([i, c, 0])

results = pd.DataFrame(data=results,
    columns=['rep', 'type', 'n'])
results.to_csv(args.out, sep='\t', index=False)
