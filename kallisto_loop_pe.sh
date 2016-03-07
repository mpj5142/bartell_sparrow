#!/bin/bash

#kallisto index -i assembly/trinity_heart.idx assembly/trinity_heart.fasta
#kallisto index -i assembly/trinity_liver.idx assembly/trinity_liver.fasta

for FILE in `ls sparrow_data/Heart*Trim_R1_P.fq.gz`
do
	FILE=`echo $FILE | perl -pe 's/\e\[?.*?[\@-~]//g'` #This line removes the auto-coloring prefixes from the file name
	FILENAME=${FILE:13}
	FILENAME=${FILENAME::10}
	echo $FILENAME
	kallisto quant -i assembly/trinity_heart.idx -o kallisto/heart/$FILENAME\_kallisto/ -b 15 -t 5 sparrow_data/$FILENAME\_R1_P.fq.gz sparrow_data/$FILENAME\_R2_P.fq.gz
done

for FILE in `ls sparrow_data/Liver*Trim_R1_P.fq.gz`
do
	FILE=`echo $FILE | perl -pe 's/\e\[?.*?[\@-~]//g'` #This line removes the auto-coloring prefixes from the file name
	FILENAME=${FILE:13}
	FILENAME=${FILENAME::10}
	echo $FILENAME
	kallisto quant -i assembly/trinity_liver.idx -o kallisto/liver/$FILENAME\_kallisto/ -b 15 -t 5 sparrow_data/$FILENAME\_R1_P.fq.gz sparrow_data/$FILENAME\_R2_P.fq.gz
done

##NOTE: Perl code for coloring prefixes came from here: 
##http://unix.stackexchange.com/questions/4527/program-that-passes-stdin-to-stdout-with-color-codes-stripped