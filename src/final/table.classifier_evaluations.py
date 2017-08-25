#!/usr/bin/env python
"""
This script writes a table with different classifier evaluation metrics.
"""
import pandas as pd
import argparse

# Add this repo to the path
import os, sys
src_dir = os.path.normpath(os.path.join(os.getcwd(), 'src/util'))
sys.path.insert(0, src_dir)
import Formatting as fmt


p = argparse.ArgumentParser()
p.add_argument('rf_in', help='path to file with tidy classifier results')
p.add_argument('fout', help='out file (Latex file)')
args = p.parse_args()

df = pd.read_csv(args.rf_in, sep='\t')
df = df.drop(['mean_fpr', 'mean_tpr'], axis=1).drop_duplicates().dropna()

# Order datasets as in Figure 1
df['dataset'] = df['dataset']\
    .replace('edd_singh', 'cdi_singh')\
    .replace('noncdi_schubert', 'cdi_schubert2')
# Need to add 'total' sample size for ordering
df['total'] = df['H_smpls'] + df['dis_smpls']
_, dataset_order = fmt.get_dataset_order(df)
df['dataset'] = df['dataset'].astype('category')
df['dataset'].cat.set_categories(dataset_order, inplace=True)
df = df.sort_values(by='dataset')

labeldict = fmt.get_labeldict(df['dataset'])
df['dataset_label'] = df['dataset'].apply(lambda x: labeldict[x])

df = df[['dataset_label', 'roc_auc', 'fisher_p', 'kappa']]
# Convert floats to reasonable strings
for col in ['roc_auc', 'fisher_p', 'kappa']:
    df[col] = df[col].apply(lambda x: '{:.2g}'.format(x))

fmt.write_latex_table(df, args.fout)
fmt.write_markdown_table(df, args.fout.rsplit('.', 1)[0] + '.md')
df.to_csv(args.fout.rsplit('.', 1)[0] + '.txt', sep='\t', index=False)
