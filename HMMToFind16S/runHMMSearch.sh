#!/bin/bash
# A simple script for scaning for scaning genomes with Hidden Markov Models.

for fasta;
do
	echo Running hmmsearch on ${fasta}
	baseFileName=$(basename ${fasta})
	hmmsearch --tblout "hmmResults$baseFileName.tsv" -A "hmmAlign$baseFileName.sto" -o "hmmHumanResults$baseFileName.txt" --acc SSU.hmm ${fasta}
done
echo All files searched.
exit 0

