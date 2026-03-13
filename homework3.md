# Homework 3
## Genome Assembly Summary
File integrity test using `md5sum`:
```bash
md5sum dmel-all-chromosome-r6.66.fasta.gz
```
Output:
ccb86e94117eb4eeaaf70efb6be1b6b9  dmel-all-chromosome-r6.66.fasta.gz

When running `faSize` on the file we get the following report:
143726002 bases (1152978 N's 142573024 real 142573024 upper 0 lower) in 1870 sequences in 1 files
Total size: mean 76858.8 sd 1382100.2 min 544 (211000022279089) max 32079331 (3R) median 1577
N count: mean 616.6 sd 6960.7
U count: mean 76242.3 sd 1379508.4
L count: mean 0.0 sd 0.0
%0.00 masked total, %0.00 masked real

1. Total number of nucleotides
```bash
    gzcat dmel-all-chromosome-r6.66.fasta.gz | grep -v ">" | wc -c
```
    There are 145523521 nucleotides.
2. Total number of Ns
```bash
gzcat dmel-all-chromosome-r6.66.fasta.gz | grep -o "N" | wc -l
```
   There are 1154851 Ns in this file.
3. Total number of sequences
   A reference fasta file is made up of 4 lines per sequence including header (chromosome + loc), sequence, and qc metrics.
```bash
gzcat dmel-all-chromosome-r6.66.fasta.gz | grep -c "^>"
```

## Summarizing Annotation file
Overview of file (first 5 lines):
```bash
bioawk -c gff '{print}' dmel-all-r6.66.gtf.gz | head -n 5
```

File integrity test using `md5sum`:
```bash
md5sum dmel-all-r6.66.gtf.gz
```
Output:
ea600dbb86f1779463f69082131753cd  dmel-all-r6.66.gtf.gz


1. Total number of features of each type, sorted from the most common to the least common
```bash
bioawk -c gff '{print $3}' dmel-all-r6.66.gtf.gz | sort | uniq -c | sort -rn
```

    Output:
    190176 exon
    163377 CDS
    46856 5UTR
    33778 3UTR
    30922 start_codon
    30862 stop_codon
    30836 mRNA
    17872 gene
    3059 ncRNA
    485 miRNA
    365 pseudogene
    312 tRNA
    270 snoRNA
    262 pre_miRNA
    115 rRNA
    32 snRNA


2. Total number of genes per chromosome arm
```bash
bioawk -c gff '$3 == "gene" {print $1}' dmel-all-r6.66.gtf.gz | sort | uniq -c | sort -rn
```

Output:
4226 3R
3649 2R
3508 2L
3481 3L
2704 X
 114 4
 113 Y
  38 mitochondrion_genome
  21 rDNA
   2 Unmapped_Scaffold_8_D1580_D1567
   2 211000022280494
   1 211000022280703
   1 211000022280481
   1 211000022280347
   1 211000022280341
   1 211000022280328
   1 211000022279681
   1 211000022279392
   1 211000022279264
   1 211000022279188
   1 211000022279165
   1 211000022278760
   1 211000022278449
   1 211000022278436
   1 211000022278279

