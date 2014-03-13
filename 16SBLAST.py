#!/usr/bin/env python 
# Created by: Lee Bergstrand
# Descript: A Bio-Python program that takes a list of query proteins and uses local BLASTp to search
#           for highly similer proteins within a local blast database (usally a local db of a target 
#           proteome). The program then BLASTps backward from the found subject protein to the proteome 
#           for which the original query protein if found in order to confirm gene orthology. 
#             
# Requirements: - This program requires the Biopython module: http://biopython.org/wiki/Download
#               - This script requires BLAST+ 2.2.9 or later.
#               - All operations are done with protien sequences.
#               - All query proteins should be from sequenced genomes in order to facilitate backwards BLAST. 
#               - MakeBlastDB must be used to create BLASTp databases for both query and subject proteomes.
#               - BLAST databases require the FASTA file they were made from to be in the same directory.
#  
# Usage: BackBLAST.py <queryGeneList.faa> <subject1.faa> 
# Example: BackBLAST.py queryGeneList.faa AUUJ00000000.faa
#----------------------------------------------------------------------------------------
#===========================================================================================================
#Imports & Setup:
import sys
import csv
import subprocess
import subprocess
from Bio import SeqIO
from Graph import Vertex
from Graph import Graph
from multiprocessing import cpu_count

processors = cpu_count() # Gets number of processor cores for BLAST.

# Dev Imports:
import time # For profiling purposes.
#===========================================================================================================
# Functions:

# 1: Checks if in proper number of arguments are passed gives instructions on proper use.
def argsCheck():
	if len(sys.argv) < 3:
		print "Orthologous Gene Finder"
		print "By Lee Bergstrand\n"
		print "Please refer to source code for documentation\n"
		print "Usage: " + sys.argv[0] + " <queryProteomes.csv> <queryGeneList.faa> <subject1.faa>\n"
		print "Examples:" + sys.argv[0] + " queryProteomes.csv queryGeneList.faa AUUJ0000000.faa"
		exit(1) # Aborts program. (exit(1) indicates that an error occured)
#-------------------------------------------------------------------------------------------------
# 2: Runs BLAST, can either be sent a fasta formatted string or a file ...
def runBLASTFor16S(query, BLASTDBFile):
	BLASTOut = subprocess.check_output(["blastn", "-db", BLASTDBFile, "-query", query, "-max_target_seqs", "1", "-num_threads", str(processors), "-outfmt", "10 qseqid sseqid length qseq evalue bitscore"]) # Runs BLASTp and save output to a string. Blastp is set to output csv which can be parsed.
	return BLASTOut
#-------------------------------------------------------------------------------------------------
# 3: Filters HSPs by Percent Identity...
def filtreBLASTCSV(BLASTOut):
	
	minIdent = 25
	
	BLASTCSVOut = BLASTOut.splitlines(True) # Converts raw BLAST csv output into list of csv rows.
	BLASTreader = csv.reader(BLASTCSVOut) # Reads BLAST csv rows as a csv.

	BLASTCSVOutFiltred = [] # Note should simply delete unwanted HSPs from curent list rather than making new list. 
					        # Rather than making a new one.
	for HSP in BLASTreader:
		if HSP[2] >= minIdent: # Filtres by minimum ident.
			# Converts each HSP parameter that should be a number to a number.
			HSP[2] = float(HSP[2]) 
			HSP[3] = float(HSP[3])
			HSP[4] = int(HSP[4])
			HSP[5] = int(HSP[5]) 
			BLASTCSVOutFiltred.append(HSP) # Appends to output array.
	
	return BLASTCSVOutFiltred
#===========================================================================================================
# Main program code:
# House keeping...
argsCheck() # Checks if the number of arguments are correct.

queryFile = sys.argv[1]

# File extension check
if not queryFile.endswith(".faa"):
	print "[Warning] " + queryFile + " may not be a amino acid fasta file!"
# File extension check
if not queryProteomesFile.endswith(".txt"):
	print "[Warning] " + queryProteomesFile + " may not be a txt file!"
	
BLASTDBFile = sys.argv[2]
print "Opening " + BLASTDBFile + "..."

print runBLASTFor16S(queryFile, BLASTDBFile)