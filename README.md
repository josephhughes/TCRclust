# TCRclust
A collection of scripts used for tracking T-cell expansions using HTS

## Primer extraction and gene usage

**SeqRenamer.pl** is a script to rename a fasta or fastq sequence either according to a text 
tab-delimited file where the user provides the column number of the old name and of the new name 
or the user provides a prefix. The prefix is used to rename the sequences along with a unique number.

**PrimerTrim.pl** is a script for counting the combinations of primers and removing the primers.

**RemoveSamePrimer.pl** is a script that does exactly what it was named as.

**FastqFilter.pl** is a script to filter the fastq reads based on minimum and maximum sequence length.

**ParseStats.pl** is a script for summarising the overall errors and the counts of forward and reverse primers. It tabulates the matrix for generating a circos plot.

## Clustering and filtering

**sort_cdhit.pl** is a script for filtering the results of cd-hit. For example, you can remove singletons or clusters that are below a certain frequency.

## Comparing clusters

**Count_clstr.pl** summarises the results of a cd-hit clustering of multiple samples by providing the unique identifiers of each sample.
