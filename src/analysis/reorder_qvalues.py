#!/usr/bin/env python
"""
Re-order rows in the given files according to order in the input tree.
Also keep only rows with at least one significant result in at least
one dataset.
"""

import argparse
import dendropy as dp
import pandas as pd

def reorder_index_from_tree(fntree, original_index):
    """
    Read a genus-level tree (tips are labeled with genus name) and reorder
    a given index according to the tree

    Parameters
    ----------
    fntree : str
        File name with newick tree. It is expected that every
        genus in the original_index is in this tree. If there is
        a genus that isn't, this code prints out the genera
        which are missing from the tree. You'll probably need to
        either update the tree with update_phyloT_wrapper.sh, or
        manually code in the missing genera in update_phyloT.py.

    original_index : list, pandas index
         List or pandas index from original df which needs to
         be re-ordered. The original_index can have OTUs with all
         taxonomic levels, or just family and genus level.

    Return
    ------
    reordered_index
    """

    # From the original index, extract just the genus name
    # and put into a dict - {genus: original_index_label}
    genus2full = {i.split(';')[-1][3:]: i for i in original_index}

    tree = dp.Tree.get(path=fntree, schema='newick')
    genera = [i.label for i in tree.taxon_namespace]
    # keep only genera from the tree that are in the original index.
    genera = [i for i in genera if i in genus2full]


    # make sure that all genera in the original_index are also in the tree
    missinggenera = [i for i in genus2full if i not in genera]
    if len(missinggenera) > 0:
        print('The following genera are missing from your tree!:')
        print('\n'.join(missinggenera))

    # Re order the original index according to the tree order
    reordered_index = [genus2full[i] for i in genera]

    return reordered_index

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
    disease_df = disease_df[['ob', 'crc', 'ibd', 'cdi', 'hiv']]
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
