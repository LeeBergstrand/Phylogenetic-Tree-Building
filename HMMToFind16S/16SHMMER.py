#!/usr/bin/env python 
# Created by: Lee Bergstrand
# Descript: A simple python program that uses HMMER to search for 16S genes within a genome. Checks
# 			both the forward and reverse strand DNA strand of the genome.          
#
# Requirements: - This program requires the Biopython module: http://biopython.org/wiki/Download
#               - This script requires HMMER 3.0 or later.
#  
# Usage: BackBLAST.py <QueryGenome.faa> <16S.hmm>
# Example: BackBLAST.py <QueryGenome.faa> AUUJ00000000.faa
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
		print "Usage: " + sys.argv[0] + "<QueryGenome.faa> <16S.hmm>\n"
		print "Examples:" + sys.argv[0] + "QueryGenome.faa <16S.hmm>"
		exit(1) # Aborts program. (exit(1) indicates that an error occured)
#-------------------------------------------------------------------------------------------------
# 2: Runs BLASTn with settings specific for extracting subject sequences.
def runHMMSearch(FASTA, HMMERDBFile):
	Found16S = True 
	process = subprocess.Popen(["hmmsearch", "--acc", "--cpu", str(processors), "-A", "tempAlign.sto", HMMERDBFile, "-"], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
	stdout, stderr = process.communicate(FASTA)
	if "No hits detected that satisfy reporting thresholds" in stdout:
		Found16S = False
	else: 
		print "16S found!"
		print stdout
	return Found16S
	
# 3: Converts the sequence record object from an RNA Stockholm file to a DNA FASTA file.
def getDNAFasta(SeqRecord):
	SeqRecord.letter_annotations = {} 	# Removes per letter annitations. Biopython throws an error if you try to 
										# reverse complement a sequence with per letter annotations (since afterward
										# these annotations would no longer be valid). We strip these per letter 
										# annotions since they are not a part of the FASTA format anyways.
	DNA = SeqRecord.seq.back_transcribe() 
	SeqRecord.seq = DNA
	return SeqRecord.format("fasta")
#===========================================================================================================
# Main program code:
# House keeping...
argsCheck(2) # Checks if the number of arguments are correct.

genomeIn = sys.argv[1]
print "Opening " + genomeIn + "..."

#fileOut = (genomeIn).strip(".fna")
#fileOut = path.split(fileOut)
#subjectAccession = fileOut[1]
#fileOut = fileOut[1] + ".16S.fna"

# File extension check
if not genomeIn.endswith(".fna"):
	print "[Warning] " + genomeIn + " may not be a nucleic acid fasta file!"
	
HMMERDBFile = sys.argv[2]
print "Opening " + HMMERDBFile + "..."
print "Searching " + genomeIn + " with " + HMMERDBFile + "..."

try:
	handle = open(genomeIn, "rU")
	SeqRecords = SeqIO.parse(handle, "fasta")
	Found16S = False
	for record in SeqRecords:
		FASTA = record.format("fasta")
		Found16S = runHMMSearch(FASTA, HMMERDBFile)
		if Found16S == True:
			handle2 = open("tempAlign.sto", "rU")
			SeqRecords = SeqIO.parse(handle2, "stockholm")
			print ">> Converting to FASTA..."
			for record in SeqRecords:
				print getDNAFasta(record)
			
	handle.close()
except IOError:
	print "Failed to open " + inFile + " or " + outFile
	exit(1)

'''
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
'''