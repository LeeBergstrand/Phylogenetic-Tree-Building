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
from Bio import SeqIO
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
#===========================================================================================================
# Main program code:
# House keeping...
argsCheck() # Checks if the number of arguments are correct.

queryFile = sys.argv[1]
print "Opening " + queryFile + "..."

# File extension check
if not queryFile.endswith(".fna"):
	print "[Warning] " + queryFile + " may not be a nucleic acid fasta file!"
	
BLASTDBFile = sys.argv[2]
print "Opening " + BLASTDBFile + "..."
print "Blasting " + queryFile + " against " + BLASTDBFile + "..."
BLASTResults = runBLASTFor16S(queryFile, BLASTDBFile)

BLASTCSVOut = BLASTResults.splitlines(True) # Converts raw BLAST csv output into list of csv rows.
BLASTreader = csv.reader(BLASTCSVOut) # Reads BLAST csv rows as a csv.

for row in BLASTreader:
	if len(row[3]) < 2000 and len(row[3]) > 1000:
		subjectAccession = row[0] 
		subject16S = row[3]
		break
		
print "Extracting 16S BLAST Results!"
FASTA = ">" + subjectAccession + " 16S rRNA\n"  + subject16S
fileOut = subjectAccession + ".16S.fna"

print "Writing results to file."
try:
	fileWriter = open(fileOut,"w")
	fileWriter.write(FASTA)
	fileWriter.close()
except IOError:
	print "Error writing " + fileOut + " to file."
print "Done."
