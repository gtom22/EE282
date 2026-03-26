#!/bin/bash 

cd /Users/gautham/Fundamentals_of_Informatics/EE282/data

# File integrity
md5sum dmel-all-chromosome-r6.66.fasta.gz

faSize dmel-all-chromosome-r6.66.fasta.gz

# get total number of nucleotides
gzcat dmel-all-chromosome-r6.66.fasta.gz | grep -v ">" | wc -c

# total number of Ns
gzcat dmel-all-chromosome-r6.66.fasta.gz | grep -o "N" | wc -l

# Total number of sequences
gzcat dmel-all-chromosome-r6.66.fasta.gz | grep -c "^>"