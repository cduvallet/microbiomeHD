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

**file-S3.core_genera.txt**

The "core" microbes that respond generally to disease. A genus is considered part of the
"core" response if it is significant (q < 0.05) in the same direction in at least two
different diseases (i.e. in at least one dataset across at least two diseases).

"Mixed" genera are those which were enriched in cases across two diseases and also
enriched in controls across two diseases.

**file-S4.literature_results.txt**

The manually curated results from each previously published paper. 
