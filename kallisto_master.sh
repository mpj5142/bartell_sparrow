#!/bin/bash

#This script creates a master Kallisto file for all transcripts across multiple samples from individual Kallisto runs.
#The output is a matrix of read counds for each transcript labeled by transcript name (col) and sample name (row).
#filemain and header are temp files, concatanated at end of script

cat Heart1Trim/abundance.tsv | sort | cut -f 1 > filemain.tsv #Get transcript names from first sample
echo "contig" > header.tsv #Start temp header file

for FILE in $(find . -name 'abundance.tsv') #Find all 'abundance' output files created by Kallisto
do
        echo $FILE #Print sample name
        echo $FILE | paste header.tsv - > temph.tsv #Get sample name and place into temp header file
        mv -f temph.tsv header.tsv
        cat $FILE | sort | cut -f 4 | paste filemain.tsv - > tempf.tsv #Get read counts and paste into temp data file
        mv -f tempf.tsv filemain.tsv

done

cat header.tsv filemain.tsv > kallisto_master_heart.tsv #Cat the temporary files into output file
rm -f header.tsv filemain.tsv #Remove temp files