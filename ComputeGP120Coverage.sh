#!/bin/bash

###########################################################################################################
# Compute average coverage of the gp120 segment for all HIV isolates used in this deep sequencing project
# To execute => ./ComputeGP120Coverage.sh in the HIV project directory which must be structured as follow;
# 
#SpecID
#   |_FASTQ
# 
# Where FASTQ for this SpecID contains his fastq.gz paired end files in format SpecID_SNNN_LNN_RN_NNN.fastq.gz.
# NNN represent any numerical value
# 
# Output is specified in the count_file variable 

count_file="/home/eric/ProjetPersonnel/ProjetsNGS/HIV/SlidingWindow/Results/Coverage/Coverage.txt"

#Average gp120 length
gp120_length=1220

#Average read length
avg_read_length=250

#Output file header
echo -e "Isolate\tAvgCoverage">>$count_file

#For all isolates fastq directories
for i in  $(ls -d *"/FASTQ/")
	do
		echo $i
		for j in $i*_R1_*.gz
                        do
                        #Extract isolate name
                        isolate=$(basename $j)
			isolate=${isolate%_*_*_*_*}
			echo "Isolate "$isolate
                        #Count number of lines in the fastq file
			count=$(zcat $j | wc -l);
                        #Divide by 4
			countf=$(($count/4));
                        #Paired end => x2
			paired_count=$(($countf*2))
                        #Compute average coverage
			cov=$((($paired_count*$avg_read_length)/$gp120_length))
			echo "Cov is "$cov
                        #Save in file
			echo -e $isolate"\t"$cov >>$count_file;
		done


done
