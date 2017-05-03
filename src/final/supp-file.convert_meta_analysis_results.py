#!/usr/bin/env python
"""
This script massages the meta-analysis results into one supplementary
file.
"""
import argparse
import pandas as pd

p = argparse.ArgumentParser()
p.add_argument('fin', help='File with genera in columns and values that '
    + ' indicate whether genera are disease (+1) or health (-1) associated '
    + ' (or mixed (0)).')
p.add_argument('fout', help='File with human-readable meta-analysis results.')
args = p.parse_args()

df = pd.read_csv(args.fin, sep='\t', index_col=0)
d = {1: 'disease', -1: 'health', 0: 'mixed'}
for val in d:
    df = df.replace(val, d[val])
df.to_csv(args.fout, sep='\t')
