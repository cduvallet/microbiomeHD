#!/bin/bash

## Given a list of genus names (or other taxonomy),
## find the NCBI ID's to then use in PhyloT

usage="get_ncbi_IDs.sh genus_list ncbi_ids [-h] -- given file with list of genus names in genus_list, fetch the NCBI IDs for use in PhyloT and save in file ncbi_ids"

if [ "$1" == "-h" ]; then
  echo $usage
  exit 0
fi

# Check if ncbi_ids file already exists.
if [ -e $2 ]; then
    echo "NCBI file ($2) exists. Over-writing."
    rm $2
fi

while read -u3 genus; do
    echo $genus
    esearch -db taxonomy -query "${genus}" | efetch -mode xml | xtract -pattern Taxon -element TaxId ScientificName Division >> $2
done 3< $1
