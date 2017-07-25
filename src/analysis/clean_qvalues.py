#!/usr/bin/env
"""
This script keeps only OTUs which were significant in at least one study.
"""
import pandas as pd
import argparse

p = argparse.ArgumentParser()
p.add_argument('qvalues_in', help='qvalues file with OTUs in rows, datasets '
               + 'in columns, and signed qvalues in values.')
p.add_argument('--qthresh', help='significant q threshold (default: '
               + ' %(default)s)', default=0.05)
args = p.parse_args()

qvalues = pd.read_csv(args.qvalues_in, sep='\t', index_col=0)
qthresh = float(args.qthresh)

# Keep only OTUs which were significant in at least one study
sig_otus = list(
           qvalues.loc[
               (qvalues.applymap(abs) < qthresh).sum(axis=1) > 0]
           .index)
qvalues = qvalues.loc[sig_otus]

# Write
newf = args.qvalues_in.rsplit('.txt', 1)[0] + '.sig.txt'
qvalues.to_csv(newf, sep='\t')
