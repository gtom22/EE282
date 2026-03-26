#!/bin/bash
  
#SBATCH --job-name=ee282       ## Name of the job.
#SBATCH -A prabhakg     ## CHANGE account to charge 
#SBATCH -p standard               ## partition name
#SBATCH --nodes=1             ## (-N) number of nodes to use
#SBATCH --ntasks=1            ## (-n) number of tasks to launch
#SBATCH --cpus-per-task=16    ## number of cores the job needs
#SBATCH --error=slurm_output/slurm-%J.err  ## error log file
#SBATCH --output=slurm_output/slurm-%J.out ## output log file

source /data/homezvol0/prabhakg/miniconda3/etc/profile.d/conda.sh
conda activate /data/homezvol0/prabhakg/miniconda3/envs/ee282


DATA=/pub/jje/ee282/ISO_HiFi_Shukla2025.fasta.gz

hifiasm -o iso1_assembly -t 16 "${DATA}"
awk '/^S/{print ">"$2"\n"$3}' iso1_assembly.bp.p_ctg.gfa > iso1_assembly.fasta

FASTA_FILE=/data/homezvol0/prabhakg/iso1_assembly.fasta

# N50 calculation
bioawk -c fastx '{print length($seq)}' "${FASTA_FILE}" | sort -rn > sequence_lengths.txt
awk '{
  len[i++] = $1;
  sum += $1
}
END {
  for (j = 0; j < i; j++){
    cumulative_sum += len[j];
    if (cumulative_sum >= sum/2){
      print len[j];
      exit
    }

  }


}' sequence_lengths.txt > n50.txt

gunzip -c dmel-all-chromosome-r6.66.fasta.gz > dmel_scaffolds.fasta
seqtk cutN -n 1 dmel_scaffolds.fasta > dmel_contigs.fasta

mkfifo pipe1 pipe2 pipe3
bioawk -c fastx '{print length($seq)}' "${FASTA_FILE}" | sort -rn > pipe1 &
bioawk -c fastx '{print length($seq)}' dmel_contigs.fasta | sort -rn > pipe2 &
bioawk -c fastx '{print length($seq)}' dmel_scaffolds.fasta | sort -rn > pipe3 &


plotCDF pipe1 pipe2 pipe3 comparison_plot.png
rm pipe1 pipe2 pipe3


compleasm run -a "${FASTA_FILE}" -l diptera -o comp_hifiasm -t 16

compleasm run -a dmel_scaffolds.fasta -l diptera -o comp_ref -t 16