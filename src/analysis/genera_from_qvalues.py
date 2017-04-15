#!/usr/bin/env python
"""
This script reads in a file with full taxonomies in the rows (and q-values
in the matrix) and writes each genus name on a line. These genus names
are fed into the esearch command to get their corresponding NCBI IDs in
downstream steps.

Full taxonomy means with the format k__;p__;...;g__
"""
import argparse
import pandas as pd

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('pvalues', help='tab-delimited file with signed '
                        + ' pvalues. Full taxonomy in index.')
    parser.add_argument('genera_outfile', help='File to write genera to or \
                        read genera from. Default is "genera.txt."', default='genera.txt')

    args = parser.parse_args()
    # Write the genera we're gonna want to get from NCBI to genera_outfile
    fnpvals = args.pvalues
    genera_outfile = args.genera_outfile

    ## Get the genera which we need in the tree and print to file
    allresults = pd.read_csv(fnpvals, sep='\t', index_col=0)
    # Keep all genera in the given dataframe
    genera = list(allresults.index)

    # Get just the genus name (i.e. no g__Genus business, just Genus)
    genera = [i.split(';')[-1][3:] for i in genera]
    with open(genera_outfile, 'w') as f:
        f.write('\n'.join(genera))
