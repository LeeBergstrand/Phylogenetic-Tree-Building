#!/bin/bash
# A simple script for searching for 16S primer sites within a genome.

for fasta;
do
	echo "Searching for primer hits in $fasta"
	echo "	Searching for the 27F Bacteria 16S rRNA Primer..."
	grep "AGAGTTTGATC[AC]TGGCTCAG" $fasta
	echo "	Searching for the 1492R Bacteria 16S rRNA Primer..."
	grep "GGTTACCTTGTTACGACTT" $fasta
	echo "	Searching for the 907R Bacteria 16S rRNA Primer..."
	grep "CCGTCAATTC[AC]TTTGAGTTT" $fasta
	echo "	Searching for the 63F Bacteria 16S rRNA Primer..."
	grep "CAGGCCTAACACATGCAAGTC" $fasta
	echo "	Searching for the 357F Bacteria 16S rRNA Primer..."
	grep "CCTACGGGAGGCAGCAG" $fasta
	echo "	Searching for the 518R Universal 16S rRNA Primer..."
	grep "ATTACCGCGGCTGCTGG" $fasta
	done
echo All files searched.
exit 0

