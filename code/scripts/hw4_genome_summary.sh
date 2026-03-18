#!/bin/bash 
# Total number of nucleotides
#Total number of Ns
#Total number of sequences

cd /Users/gautham/Fundamentals_of_Informatics/EE282/data

# partition the genome into 100kb and smaller and larger than 100kb
faSize -detailed dmel-all-chromosome-r6.66.fasta.gz | bioawk -v OFS="\t" '{if ($2 < 100000) print $0, "smaller"; else print $0, "larger"}' > dmel-all-chromosome-r6.66.size_summary.txt
# create separate files for each partition
awk -F "\t" '$3 == "smaller" {print $1}' dmel-all-chromosome-r6.66.size_summary.txt > smaller_ids.txt
awk -F "\t" '$3 == "larger" {print $1}' dmel-all-chromosome-r6.66.size_summary.txt > larger_ids.txt


seqtk subseq dmel-all-chromosome-r6.66.fasta.gz smaller_ids.txt | gzip > dmel.smaller.fasta.gz
seqtk subseq dmel-all-chromosome-r6.66.fasta.gz larger_ids.txt | gzip > dmel.larger.fasta.gz

for file in dmel.smaller.fasta.gz dmel.larger.fasta.gz
do
    base=$(basename "$file" .fasta.gz)  
    faSize "$file"
    
    bioawk -c fastx -v name="$base" '{l = length($seq); 
        g = gc($seq);
        n = gsub(/[Nn]/, "", $seq); 
        gc_bases = (l * g);
        print $name, l, n, g;
        total_l += l; total_n += n; 
        total_gc_bases += gc_bases} 
        END  {
            print "Summary" > "/dev/stderr";
            print "Total bases: " total_l > "/dev/stderr";
            print "Total Ns: " total_n > "/dev/stderr";
            print "Overall GC: " (total_gc_bases/total_l) > "/dev/stderr";
        }' "$file" > "${base}_stats.txt"


    cut -f2  "${base}_stats.txt" | sort -rn > "${base}_lengths.txt"
    plotCDF "${base}_lengths.txt" "${base}_plot.png"
    Rscript ../code/scripts/plothisto.R "${base}_stats.txt" "$base"
done



