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
p.add_argument('qvalues_in', help='qvalues file with OTUs in rows, datasets '
               + 'in columns, and signed qvalues in values.')
p.add_argument('meta_sig', help='file with disease-wise significant bugs')
p.add_argument('overall_sig', help='file with overall significant bugs')
p.add_argument('fntree', help='path to tree file (newick)')
p.add_argument('--qthresh', help='significant q threshold (default: '
               + ' %(default)s)', default=0.05)

args = p.parse_args()

qvalues = pd.read_csv(args.qvalues_in, sep='\t', index_col=0)
disease_df = pd.read_csv(args.meta_sig, sep='\t', index_col=0)
overall_df = pd.read_csv(args.overall_sig, sep='\t', index_col=0)
qthresh = float(args.qthresh)

sig_otus = list(\
           qvalues.loc[\
               (qvalues.applymap(abs) < qthresh).sum(axis=1) > 0]\
           .index)

# Re-order these significant OTUs according to the tree
sig_otus = reorder_index_from_tree(args.fntree, sig_otus)

## Update the other meta-analysis dataframes to have the same rows
overall_df = overall_df.loc[sig_otus]
qvalues = qvalues.loc[sig_otus]
disease_df = disease_df.loc[sig_otus]

newf = args.qvalues_in.rsplit('.txt', 1)[0] + '.sig_ordered.txt'
qvalues.to_csv(newf, sep='\t')
newf = args.meta_sig.rsplit('.txt', 1)[0] + '.sig_ordered.txt'
disease_df.to_csv(newf, sep='\t')
newf = args.overall_sig.rsplit('.txt', 1)[0] + '.sig_ordered.txt'
overall_df.to_csv(newf, sep='\t')
