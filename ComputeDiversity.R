###############################################
# This script compute nucleotidic diversity   #
# from an alignment in fasta format           #
# Execute it as follow;                       #
# Rscript ComputeDiversity.R <alignment.fas>  #
###############################################

#read the fasta alignment
args=commandArgs(TRUE)
my_fasta=args[1]

#For diversity calculation
library(pegas)

#For sequence alignment
library(msa)

#For biological sequence manipulation
library(seqinr)

#Non-aligned fasta file
non_align=readDNAStringSet(my_fasta)

#Do alignment
align=msa(non_align,"ClustalW",type="dna")
#Convert to DNAbin class prior to compute diversity
align_bin=as.DNAbin(align)

#Compute diversity
nuc.div(align_bin)

