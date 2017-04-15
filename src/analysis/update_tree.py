#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Script to update the phyloT tree with some manually coded genera
which did not have an NCBI ID using esearch.

@author: claire
"""
import pandas as pd
import dendropy as dp
import argparse

def hard_coded_additions():
    """
    Hard coded dictionary with {nearest_parent: [list, of, missing, genera]}

    For example, any genus belonging to the family Clostridiaceae would be in:
        {'Clostridiaceae': ['Clostridium sensu stricto', 'Clostridium_XI', ...]}

    """

    genera_dict = {'Clostridiaceae':  ["Clostridium_sensu_stricto",
                                       'Anaerobacter'],
                   'Lachnospiraceae': ['Clostridium_XlVa',
                                       'Clostridium_XlVb',
                                       'Lachnospiracea_incertae_sedis',
                                       'Ruminococcus2'],
                    'Ruminococcaceae': ['Clostridium_III',
                                        'Clostridium_IV'],
                    'Erysipelotrichaceae': ['Clostridium_XVIII',
                                            'Erysipelotrichaceae_incertae_sedis'],
                    'Peptostreptococcaceae' : ['Clostridium_XI'],
                    'Enterobacteriaceae': ['Escherichia/Shigella'],
                    'Prevotellaceae': ['Xylanibacter'],
                    'Flavobacteriaceae': ['Planobacterium'],
                    'Fusobacteriaceae' : ['Clostridium_XIX'],
                    'Clostridiales': ['Clostridium_XII'],
                    'Bacteria': ["Saccharibacteria_genera_incertae_sedis",
                                 "Subdivision5_genera_incertae_sedis",
                                 'GpI']}

    return genera_dict

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('genera_file', help='File with genera on each line')
    parser.add_argument('tree_in', help='Input phyloT tree (newick)')
    parser.add_argument('tree_out', help='Output tree, with any genera '
                        + 'in the given genera_file that are not in the '
                        + 'given tree added manually by hard-coded dict '
                        + 'in this script.')

    args = parser.parse_args()

    # Read in genera in the data file
    data_genera = [i.strip() for i in open(args.genera_file, 'r').readlines()]

    # Read in the tree
    tree = dp.Tree.get(path=args.tree_in, schema='newick')

    # Find the genera which should be in the tree but which aren't
    tree_genera = [i.label for i in tree.taxon_namespace]
    missing_genera = [i for i in data_genera if i not in tree_genera]

    # Get the hard-coded additions (like, where each missing genus fits)
    genera_dict = hard_coded_additions()

    # Check that I have all the missing_genera in my genera_dict
    added = [i for sublist in genera_dict.values() for i in sublist]
    toadd = [i for i in missing_genera if i not in added]
    if len(toadd) > 0:
        print('Warning: you are missing the following genera from your manual update!!')
        print('\n'.join(toadd))
        raise ValueError('Missing genera will not go quietly.')

    # Add the additions (see code below woops)
    for parent in genera_dict:
        for child in genera_dict[parent]:
            parent_node = tree.find_node_with_label(parent)
            #print(parent_node)
            newchild = dp.Taxon(label=str(child))
            parent_node.new_child(taxon=newchild)
    tree.reconstruct_taxon_namespace()
    tree.write(path=args.tree_out, schema='newick')
    # FYI useful dendropy resources:
    #https://groups.google.com/forum/#!topic/dendropy-users/Xy5KxHFur8E
    #http://stackoverflow.com/questions/26835804/adding-new-nodes-to-a-tree-by-dendropy
