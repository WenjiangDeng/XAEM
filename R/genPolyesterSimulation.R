# Example of generating simulated data of transcripts using polyester packages
rm(list=ls())

args = commandArgs(trailingOnly=TRUE)
cat("\nNumber of arguments: ",length(args))
cat("\nList of arguments: ",args,"\n")

fasta_file=as.character(args[1])
outdir = as.character(args[2])


## try http:// if https:// URLs are not supported
source("http://bioconductor.org/biocLite.R")
biocLite("polyester")
source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings")

library(polyester)
library(Biostrings)
fasta = readDNAStringSet(fasta_file)
# generate reads with coverage of 2
unif.countmat=2*width(fasta)
# generate at least 1000 reads
unif.countmat[unif.countmat<1000]=1000
unif.countmat=as.matrix(unif.countmat)
simulate_experiment_countmat(fasta_file, readmat=unif.countmat, outdir=outdir,error_rate=0.0, strand_specific=FALSE) 
rm(list=ls())

