#!/usr/bin/env python 
# Created by: Lee Bergstrand
# Descript: A simple python program that use BLASTn to search for 16S genes within a geneome. 
#			Extracts the aligned query sequence from the BLAST results (This should be the 16S gene in the query genome).              
#
# Requirements: - This program requires the Biopython module: http://biopython.org/wiki/Download
#               - This script requires BLAST+ 2.2.9 or later.
#               - MakeNABlastDB must be used to create BLASTn databases for both query and subject proteomes.
#               - BLAST databases require that the FASTA file they were made from be in the same directory.
#  
# Usage: BackBLAST.py <queryGeneList.faa> <subject1.faa> 
# Example: BackBLAST.py queryGeneList.faa AUUJ00000000.faa
#----------------------------------------------------------------------------------------
#===========================================================================================================
#Imports & Setup:
import sys
import csv
from os import path 
import subprocess
from Bio import SeqIO
from multiprocessing import cpu_count

processors = cpu_count() # Gets number of processor cores for BLAST.

# Dev Imports:
import time # For profiling purposes.
#===========================================================================================================
# Functions:

# 1: Checks if in proper number of arguments are passed gives instructions on proper use.
def argsCheck(argsCount):
	if len(sys.argv) < 3:
		print "Orthologous Gene Finder"
		print "By Lee Bergstrand\n"
		print "Please refer to source code for documentation\n"
		print "Usage: " + sys.argv[0] + "<QueryGenome.faa> <16SDataBase.faa>\n"
		print "Examples:" + sys.argv[0] + "QueryGenome.faa .faa"
		exit(1) # Aborts program. (exit(1) indicates that an error occured)
#-------------------------------------------------------------------------------------------------
# 2: Runs BLASTn with settings specific for extracting subject sequences.
def runBLASTFor16S(query, BLASTDBFile):
	# Runs BLASTn and saves the output to a string. Blastn is set to output a csv which can be parsed by Pythons CSV module.
	BLASTOut = subprocess.check_output(["blastn", "-db", BLASTDBFile, "-query", query, "-max_target_seqs", "1", "-num_threads", str(processors), "-outfmt", "10 qseqid sseqid length qseq evalue bitscore"]) 
	return BLASTOut
#===========================================================================================================
# Main program code:
# House keeping...
argsCheck(2) # Checks if the number of arguments are correct.

queryFile = sys.argv[1]
print "Opening " + queryFile + "..."

fileOut = (queryFile).strip(".fna")
fileOut = path.split(fileOut)
subjectAccession = fileOut[1]
fileOut = fileOut[1] + ".16S.fna"

# File extension check
if not queryFile.endswith(".fna"):
	print "[Warning] " + queryFile + " may not be a nucleic acid fasta file!"
	
BLASTDBFile = sys.argv[2]
print "Opening " + BLASTDBFile + "..."
print "Blasting " + queryFile + " against " + BLASTDBFile + "..."
BLASTResults = runBLASTFor16S(queryFile, BLASTDBFile)

if not BLASTResults: # If there are no BLAST results aborts the program.
	print "\nUnfortunetly there are no 16S BLAST results for " + queryFile + ". Try using another BLAST DB. or"
	print "you may also want to try another method to find 16S other than BLAST (eg. HMMs).\n"
	exit(1)  # Aborts program. (exit(1) indicates that an error occured)
	
BLASTCSVOut = BLASTResults.splitlines(True) # Converts raw BLAST csv output into list of csv rows.
BLASTreader = csv.reader(BLASTCSVOut) # Reads BLAST csv rows as a csv.

Found16S = False
for row in BLASTreader:
	if len(row[3]) < 2000 and len(row[3]) > 1000: # 16S genes are around 1500 B.P. This fitres out partial sequence or really large sequences.
		subject16S = row[3]
		Found16S = True
		break

if Found16S == False: # If there are no BLAST that are around the size of 16S rRNA abort the program.
	print "\nUnfortunetly there are no 16S BLAST results for " + queryFile + ". Try using another BLAST DB or"
	print "you may also want to try another method to find 16S other than BLAST (eg. HMMs).\n"
	exit(1)  # Aborts program. (exit(1) indicates that an error occured)

print "Extracting 16S BLAST Results!"
FASTA = ">" + subjectAccession + " 16S rRNA\n"  + subject16S

print "Writing results to file."
try:
	fileWriter = open(fileOut,"w")
	fileWriter.write(FASTA)
	fileWriter.close()
except IOError:
	print "Error writing " + fileOut + " to file."
print "Done.\n"
