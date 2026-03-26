#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) > 1) {
  stop("\nUsage: plothisto infile(s) outfile\nTakes one or more contig length files in descending order\nand returns a single png figure of the histogram of them.")
}


infile <- args[1]
prefix <- args[2]

data <- read.table(infile, header = FALSE)
colnames(data) <- c("name", "length", "ns", "gc_frac")

png(filename = paste0(prefix, "_length_hist.png"), width = 640, height = 480)
hist(log10(data$length), 
     breaks = 50, 
     col = "blue", 
     border = "white",
     main = paste("Sequence Length Distribution:", prefix),
     xlab = "Log10 Sequence Length",
     ylab = "Frequency")
dev.off()


png(filename = paste0(prefix, "_gc_hist.png"), width = 640, height = 480)
hist(data$gc_frac * 100, 
     breaks = 50, 
     col = "red", 
     border = "white",
     main = paste("GC Content Distribution:", prefix),
     xlab = "GC Percentage",
     ylab = "Frequency")
dev.off()
