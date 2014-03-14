#!/bin/bash
# A simple script for the batch creation of blast databases from a directory with fasta files.

for fasta;
do
	echo Running blast on $fasta
	python 16SBLAST.py $fasta ./Example16DB/RDPActinoBacteria16S.fna
	done
echo All files BLASTed.
exit 0

