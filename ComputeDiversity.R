args=commandArgs(TRUE)
my_fasta=args[1]

library(pegas)
library(msa)
library(seqinr)

non_align=readDNAStringSet(my_fasta)
align=msa(non_align,"ClustalW",type="dna")
align_bin=as.DNAbin(align)
nuc.div(align_bin)

