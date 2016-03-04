#!/bin/bash

makeblastdb -in transps/zebra_finch_protein.fa -dbtype prot

blastx -db transps/zebra_finch_protein.fa -query assembly/trinity_heart.fasta -out transps/blastx_heart.out -outfmt 6 -evalue 0.01 -max_target_seqs 20

perl /home/matt/work/src/TransPS1.1.0/transps.pl -t assembly/trinity_heart.fasta -b transps/blastx_heart.out