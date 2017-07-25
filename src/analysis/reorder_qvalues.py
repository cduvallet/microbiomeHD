#!/usr/bin/env python
"""
Re-order rows in the given files according to order in the input tree.
Also keep only rows with at least one significant result in at least
one dataset.
"""

import argparse
import pandas as pd

# Add src/util/ to the path
import os, sys
src_dir = os.path.normpath(os.path.join(os.getcwd(), 'src/util'))
sys.path.insert(0, src_dir)
from Formatting import reorder_index_from_tree

p = argparse.ArgumentParser()
p.add_argument('--disease-df', help='file with disease-wise significant bugs')
p.add_argument('--overall', help='file with overall significant bugs')
p.add_argument('--qvalues', help='file with clean qvalues file. If this is '
    + 'given, the input file(s) will keep only the rows also in this file.')
p.add_argument('--do-qvals', help='whether to reorder the given qvalues file.',
    action='store_true')
p.add_argument('fntree', help='path to tree file (newick)')
args = p.parse_args()

if args.qvalues is not None:
    sig_otus = pd.read_csv(args.qvalues, sep='\t', index_col=0).index
elif args.disease_df is not None:
    sig_otus = pd.read_csv(args.disease_df, sep='\t', index_col=0).index
elif args.overall is not None:
    sig_otus = pd.read_csv(args.overall, sep='\t', index_col=0).index

# Re-order rows according to the tree
sig_otus = reorder_index_from_tree(args.fntree, sig_otus)

if args.do_qvals and args.qvalues:
    qvalues = pd.read_csv(args.qvalues, sep='\t', index_col=0)
    qvalues = qvalues.loc[sig_otus]
    newf = args.qvalues.rsplit('.txt', 1)[0] + '_ordered.txt'
    qvalues.to_csv(newf, sep='\t')

if args.disease_df:
    disease_df = pd.read_csv(args.disease_df, sep='\t', index_col=0)
    # Re-order rows
    disease_df = disease_df.loc[sig_otus]
    # And re-order columns manually
    disease_df = disease_df[['cdi', 'ob', 'crc', 'ibd', 'hiv']]
    # Write to file
    newf = args.disease_df.rsplit('.txt', 1)[0] + '.sig_ordered.txt'
    disease_df.to_csv(newf, sep='\t')

if args.overall:
    overall_df = pd.read_csv(args.overall, sep='\t', index_col=0)
    # Reorder rows
    overall_df = overall_df.loc[sig_otus]
    # Write to file
    newf = args.overall.rsplit('.txt', 1)[0] + '.sig_ordered.txt'
    overall_df.to_csv(newf, sep='\t')
