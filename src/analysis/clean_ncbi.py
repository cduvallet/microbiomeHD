#!/usr/bin/env python
"""
This script takes in the output of the esearch command (a tab-delimited
file with id, genus, and kingdom names in columns) and writes a file
with just the Bacterial genus IDs.

Some of the genera return non-Bacterial things when querying against
NCBI, so that's why this step is necessary.
"""
import os
import argparse
import pandas as pd

p = argparse.ArgumentParser()
p.add_argument('ncbi_in', help='NCBI file with id, genus, and kingdom columns')
p.add_argument('ncbi_clean', help='out file with only Bacteria NCBI entries')
p.add_argument('ncbi_out', help='out file with just NCBI IDs on each line')
args = p.parse_args()


## Clean up the NCBI IDs
ids = pd.read_csv(args.ncbi_in, sep='\t', names=['id', 'genus', 'kingdom'])

# Drop duplicates
ids = ids.drop_duplicates()

## Find returned genera which aren't bacterial
not_bacteria = ids.query('kingdom != "Bacteria"')

drop_rows = []
for g, subdf in not_bacteria.groupby('genus'):
    # Then there are two entries with that genus name, see if there's one
    # that's bacteria
    tmp = ids.query('genus == @g').query('kingdom == "Bacteria"')
    #print(tmp)
    if tmp.shape[0] > 0:
        drop_rows += subdf.index.tolist()
    else:
        print(subdf)
# Drop rows corresponding to genera which do have a Bacterial ID
# (and which have duplicate in other kingdoms)
ids = ids.drop(drop_rows)

# Overwrite NCBI file
ids.to_csv(args.ncbi_clean, sep='\t', index=False)

## Print just the NCBI IDs to a file for input into phyloT
# grab all but the first line | grab first column > save to file
cmd = "tail -n +2 {} | cut -f 1 > {}"\
            .format(args.ncbi_clean, args.ncbi_out)
os.system(cmd)
# Add Akkermansia, for some reason my query script doesn't find it
cmd = 'echo "239934" >> {}'.format(args.ncbi_out)
os.system(cmd)
