#!/usr/bin/env python
"""
This script lists the significant bugs for a given dataset in tex format.
"""

import argparse
import pandas as pd
import numpy as np

p = argparse.ArgumentParser()
p.add_argument('qvalues', help='file with signed qvalues')
p.add_argument('dataset', help='dataset ID, as in qvalues file columns')
p.add_argument('--qthresh', help='significance threshold [default: '
    + '%(default)s]', type=float, default=0.05)
args = p.parse_args()

qvals = pd.read_csv(args.qvalues, sep='\t', index_col=0)
qthresh = args.qthresh

dirs = {-1: 'healthy', 1: 'disease'}
sigmap = lambda x: dirs[np.sign(x)] if abs(x) < qthresh and abs(x) > 0 \
    else np.nan

sig = qvals[args.dataset].apply(sigmap)

healthy = [i.split(';')[-1][3:] for i in
    sig[sig.astype(str) == 'healthy'].index]
disease = [i.split(';')[-1][3:] for i in
    sig[sig.astype(str) == 'disease'].index]

print('\nhealthy')
print(', '.join(['\\textit{' + '\_'.join(i.split('_')) + '}' for i in healthy]))
print('\ndisease')
print(', '.join(['\\textit{' + '\_'.join(i.split('_')) + '}' for i in disease]))
