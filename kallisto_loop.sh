#!/bin/bash

kallisto index -i assembly/trinity_heart.idx assembly/trinity_heart.fasta
kallisto index -i assembly/trinity_liver.idx assembly/trinity_liver.fasta

for FILE in `ls SE_Reads/concat/*Heart.fastq`
do
	FILENAME=${FILE:16}
	kallisto quant -i assembly/trinity_heart.idx -o kallisto/$FILENAME\_kallisto/ --single -l 200 -s 1 -b 15 -t 5 $FILE
done

for FILE in `ls SE_Reads/concat/*Liver.fastq`
do
	FILENAME=${FILE:16}
	kallisto quant -i assembly/trinity_liver.idx -o kallisto/$FILENAME\_kallisto/ --single -l 200 -s 1 -b 15 -t 5 $FILE
done