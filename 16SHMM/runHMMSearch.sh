#!/bin/bash
# A simple script for the batch creation of blast databases from a directory with fasta files.

for fasta;
do
	echo Running hmmsearch on $fasta
	baseFileName=$(basename $fasta)
	hmmsearch --tblout "hmmResults$baseFileName.tsv" -A "hmmAlign$baseFileName.stk" -o "hmmHumanResults$baseFileName.txt" --acc SSU.hmm $fasta
done
echo All files searched.
exit 0

