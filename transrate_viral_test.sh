#Test of the TransRate transcriptome contig scoring software
#Commands run by M. Jensen, Dr. Istvan Albert's lab, 2/9/2016

#For this test, I downloaded a small viral genome (West Nile virus), simulated a transcriptome read set #from the genome, and ran it through a de novo transcriptome pipeline. 

#Download West Nile Virus genome from NCBI
efetch -db=nuccore -format=fasta -id=AJ965628 > AJ965628.fa

#Simulate Illumina 10X coverage of the genome.
# N= 10977 bp/100bp ~ 110 reads *10X coverage = 1100 simulated reads
# Generate read sets error-free (as control), and with 1%, 5%, 10%, and 20% error.

wgsim -N 1100 -e 0.0 -1 100 -2 100 AJ965628.fa sim_error0_1.fq sim_error0_2.fq
wgsim -N 1100 -e 0.01 -1 100 -2 100 AJ965628.fa sim_error1_1.fq sim_error1_2.fq
wgsim -N 1100 -e 0.05 -1 100 -2 100 AJ965628.fa sim_error5_1.fq sim_error5_2.fq
wgsim -N 1100 -e 0.1 -1 100 -2 100 AJ965628.fa sim_error10_1.fq sim_error10_2.fq
wgsim -N 1100 -e 0.2 -1 100 -2 100 AJ965628.fa sim_error20_1.fq sim_error20_2.fq

#Assemble reads with Velvet
#First use velveth to construct assembly datasets for velvetg, and then use velvetg to generate assembly

velveth sim_error_0/ 21 -shortPaired -fastq sim_error0_1.fq sim_error0_2.fq
velvetg sim_error_0/ -ins_length 250 -exp_cov 10.0
velveth sim_error_1/ 21 -shortPaired -fastq sim_error1_1.fq sim_error1_2.fq
velvetg sim_error_1/ -ins_length 250 -exp_cov 10.0
velveth sim_error_5/ 21 -shortPaired -fastq sim_error5_1.fq sim_error5_2.fq
velvetg sim_error_5/ -ins_length 250 -exp_cov 10.0
velveth sim_error_10/ 21 -shortPaired -fastq sim_error10_1.fq sim_error10_2.fq
velvetg sim_error_10/ -ins_length 250 -exp_cov 10.0
velveth sim_error_20/ 21 -shortPaired -fastq sim_error20_1.fq sim_error20_2.fq
velvetg sim_error_20/ -ins_length 250 -exp_cov 10.0

#To get stats on each of the resulting "contigs.fa" files...
cat contigs.fa | grep ">" | wc -l #Get the total number of contigs
cat contigs.fa | grep ">" | cut -d '_' -f 4 | sort -nr | head -1 #Get size of largest contig
cat contigs.fa | grep ">" | cut -d '_' -f 4 | paste -sd+ - | bc #Get the sum of all contig lengths

#DATA
# Error Rate 	# Contigs 	Longest Contig 		Total length
#	0%			1 			10946 				10946
#	1%			56			3193				11387
#	5%			596			268					27284
#	10%			941			244					50307
#	20%			92	 		946 				6332

#Run Transrate on each assembly independently
#(Note: .fastq files were moved into Velvet output directories before this step)
transrate --assembly contigs.fa --left sim_error0_1.fq --right sim_error0_2.fq 
transrate --assembly contigs.fa --left sim_error1_1.fq --right sim_error1_2.fq 
transrate --assembly contigs.fa --left sim_error5_1.fq --right sim_error5_2.fq
transrate --assembly contigs.fa --left sim_error10_1.fq --right sim_error10_2.fq 
transrate --assembly contigs.fa --left sim_error20_1.fq --right sim_error20_2.fq 

