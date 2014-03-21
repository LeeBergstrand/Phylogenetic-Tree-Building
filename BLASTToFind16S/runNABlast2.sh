#!/bin/bash
# A simple script for the batch creation of blast databases from a directory with fasta files.

for fasta;
do
	echo Running blast on $fasta
	blastn -db ./Example16DB/RDPActinoBacteria16S.fna -query $fasta -outfmt "10 qseqid sseqid length qseq evalue bitscore"
done
echo All files BLASTed.
exit 0

