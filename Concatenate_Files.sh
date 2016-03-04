 #Script to concatenate all reads from a tissue to a single fastq file (ONLY DO THIS FOR FOLLOWUP, SINGLE-END READS!)
#WJH 01/29/16
#Modified by MCJ 03/02/16

#Directories

BASE_DIR=~/work/data/sparrow/SE_Reads/ #Change to server location...

RUN1=160121_7001126F_0084_AHGVYMBCXX/fastq
RUN2=160122_7001126F_0085_AHGTLCBCXX/fastq
RUN3=160122_7001126F_0086_BHH35VBCXX/fastq
RUN4=160125_7001126F_0087_AHH2Y5BCXX/fastq
RUN5=160125_7001126F_0088_BHH2WVBCXX/fastq
RUN6=160126_7001126F_0089_AHH552BCXX/fastq

OUTPUT_DIR=~/work/data/sparrow/SE_Reads/concat #Change to server location...

#Loop to find all files and then concatenate them all together from all 6 separate runs:
cd $BASE_DIR/$RUN1

NAMES=$(find -name *.fastq.gz)

for n in $NAMES
	do 
	FULL=${n#./}
	PREFIX=${n%_R1.fastq*}
	PREFIX=${PREFIX##*Sample_Sparrow}
	zcat $BASE_DIR/$RUN1/$FULL $BASE_DIR/$RUN2/$FULL $BASE_DIR/$RUN3/$FULL $BASE_DIR/$RUN4/$FULL $BASE_DIR/$RUN5/$FULL $BASE_DIR/$RUN6/$FULL > $OUTPUT_DIR/${PREFIX}.fastq
done