# Lapels #

Lapels is used to remap in silico (pseudo) genome alignments, in the form of a BAM file, back to the reference sequence. This entails the removal of all indels (via the cigar string modifications, the underlying sequence is unaltered) and adjustment fragment and mate starting positions. Lapels also annotates the number and types (SNPs, insertions, and deletions) of allele variants seen in each read.

The input includes the BAM file of in silico (pseudo) genome alignments and the MOD file associated with the FASTA sequences used in the alignment.

Specifically, the MOD and FASTA files of 8 CC founders, 12 CC F1 strains, as well as strains from Wellcome Trust Sanger Institute, can be downloaded from:
http://csbio.unc.edu/CCstatus/index.py?run=Pseudo .(Please bundle MOD and FASTA while downloading.)

The output is a BAM file with corrected reads positions, adjusted cigar strings, and annotated tags. It has been tested to be compatible with downstream tools, such as [IGV](http://www.broadinstitute.org/igv/) (using the reference genome) and [Cufflinks](http://cufflinks.cbcb.umd.edu/) (using any referenced based transcript library).

The code for Lapels is written in Python. It requires the [pysam](http://code.google.com/p/pysam/) library and the [argparse](http://code.google.com/p/argparse/) library.


https://pypi.python.org/pypi/lapels