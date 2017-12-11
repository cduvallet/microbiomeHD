# Supplementary Files

**file-S1.qvalues.txt**

This file contains the qvalues for case-control comparisons for each genus in each dataset.
Negative values indicate that the mean abundance in controls was higher than in cases.
Positive values indicate that the mean abundance in cases was higher that in controls.
Missing values mean that that genus was not present in that dataset.

**file-S2.disease_specific_genera.txt**

Genera which are consistently associated with individual diseases, for the diseases which
had at least 3 independent datasets.

Genera are considered disease-associated if they are significant (q < 0.05) in the same
direction in at least two different datasets of the same disease.

Note that if a genus is significantly health-associated in 2 datasets and signficantly
disease-associated in 1 dataset of the same disease, then it is *not* considered to be
consistently associated with that disease. In other words, being both health- and
disease-associated cancels out and only the *net* association is considered.

**file-S3.nonspecific_genera.txt**

The microbes that respond non-specifically to disease. A genus is considered part of the
shared response if it is significant (q < 0.05) in the same direction in at least two
different diseases (i.e. in at least one dataset across at least two diseases).

"Mixed" genera are those which were enriched in cases across two diseases and also
enriched in controls across two diseases.

**file-S4.literature_results.txt**

The manually curated results from each previously published paper.

**file-S5.effects.txt**

Logfold difference of mean abundance in cases and controls.
For each genus, this is calculated as:
`log2(mean_abundance_in_cases/mean_abundance_in_controls)`.
If the mean abundance in both cases and controls is 0, this value
is 0. If the mean abundance in controls is 0 and in cases is not 0,
the maximum value of the table is filled in. If the mean abundance
in cases is 0 and in controls is not 0, the minimum value of the
table is filled in.
