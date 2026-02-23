#!/bin/bash 

cd /Users/gautham/Fundamentals_of_Informatics/EE282/data

# overview
bioawk -c gff '{print}' dmel-all-r6.66.gtf.gz | head -n 5

# File integrity
md5sum dmel-all-r6.66.gtf.gz

# Total number of features of each type, sorted from the most common to the least common
bioawk -c gff '{print $3}' dmel-all-r6.66.gtf.gz | sort | uniq -c | sort -rn

# Total number of genes per chromosome arm
bioawk -c gff '$3 == "gene" {print $1}' dmel-all-r6.66.gtf.gz | sort | uniq -c | sort -rn