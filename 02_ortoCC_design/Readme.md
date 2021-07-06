# Designing CCNG1

This directory contains the code used to design the orthogonal sets for the CCNG1 orthogonal sets. In the code and old nodes this library is also  referred to as R8000, because it has about 8000 pairs.

The notebooks contain the code and comments need to run. They are numbered in the order of execution (with at two digit prefix).

The order is:

* 10-all.ipynb  -- design with  all posible AA positions
* 20-1or2-N.ipynb -- allow a maximum of two ASP
* 21-1or2-N-not_first.ipynb -- allow a maximum of two ASN but not in the first heptade
* 23-1or2N-all-same_GE.ipynb -- have same amino acid in the g and e position of the same heptade
* 30-1or2-N-not_first_only.ipynb -- allow a maximum of two ASN but not in the first heptade if there is only one ASN
* 31-1or2-N-not_first_only-same_GE.ipynb -- allow a maximum of two ASN but not in the first heptade if there is only one ASN and keep the same ge residues
* 40-score-all.ipynb -- re-scores all the sets
* 50-check-sets.ipynb -- checks set statistics and creates nice excel tables for each set
* 60-make_bcf.ipynb -- makes sequences with different bcf positions (different backgorunds)
* 70-collect_pairs.ipynb -- collect the sequences of all sets into a fasta file
