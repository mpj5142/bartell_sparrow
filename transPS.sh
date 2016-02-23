#!/bin/bash

makeblastdb -in transps/white_sparrow_prot.fa -dbtype prot

blastx -db transps/white_sparrow_prot.fa -query sparrow_data/trinity_output_heart/Trinity.fasta -out transps/blastx_heart.out -outfmt 6 -evalue 0.01 -max_target_seqs 20

perl /home/matt/work/src/TransPS1.1.0/transps.pl -t sparrow_data/trinity_output_heart/Trinity.fasta -b transps/blastx.out