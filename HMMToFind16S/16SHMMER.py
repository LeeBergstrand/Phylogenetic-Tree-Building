#!/usr/bin/env python 
# Created by: Lee Bergstrand
# Descript: A simple python program that uses HMMER to search for 16S genes within a genome. Checks
# 			both the forward and reverse strand DNA strand of the genome.          
#
# Requirements: - This program requires the Biopython module: http://biopython.org/wiki/Download
#               - This script requires HMMER 3.0 or later.
#  
# Usage: BackBLAST.py <Querygenome.faa> <16S.hmm>
# Example: BackBLAST.py <Querygenome.faa> AUUJ00000000.faa
#----------------------------------------------------------------------------------------
#===========================================================================================================
#Imports & Setup:
import sys
import csv
import cStringIO
from os import path 
import subprocess
from Bio import SeqIO
from multiprocessing import cpu_count

processors = cpu_count() # Gets number of processor cores for HMMER.

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
		print "Usage: " + sys.argv[0] + "<Querygenome.faa> <16S.hmm>\n"
		print "Examples:" + sys.argv[0] + "Querygenome.faa <16S.hmm>"
		exit(1) # Aborts program. (exit(1) indicates that an error occured)
#-------------------------------------------------------------------------------------------------
# 2: Runs BLASTn with settings specific for extracting subject sequences.
def runHMMSearch(FASTA, HMMERDBFile):
	Found16S = True
	process = subprocess.Popen(["hmmsearch", "--acc", "--cpu", str(processors), "-A", "tempAlign.sto", HMMERDBFile, "-"], stdin = subprocess.PIPE, stdout = subprocess.PIPE, bufsize = 1)
	stdout = process.communicate(FASTA)[0] # This returns a list with both stderr and stdout. Only want stdout whic is first element.
	if "No hits detected that satisfy reporting thresholds" in stdout:
		Found16S = False
	return Found16S
#-------------------------------------------------------------------------------------------------
# 3: Converts the sequence record object from an RNA Stockholm file to a DNA FASTA file.
def getDNAFasta(SeqRecord):
	SeqRecord.letter_annotations = {} 	# Removes per letter annitations. Biopython throws an error if you try to 
										# reverse complement a sequence with per letter annotations (since afterward
										# these annotations would no longer be valid). We strip these per letter 
										# annotions since they are not a part of the FASTA format anyways.
	DNA = SeqRecord.seq.back_transcribe() 
	SeqRecord.seq = DNA
	FASTA = SeqRecord.format("fasta") 
	FASTA = fastaClean(FASTA)
	return FASTA
#-------------------------------------------------------------------------------------------------
# 4: Addes sequence files to the lists.
def add16SSequences(SixteenSSubunits):
	handle2 = open("tempAlign.sto", "rU") # Read in the alignment file created by runHMMSearch.
	SixTeens = SeqIO.parse(handle2, "stockholm") # Parse the alignment file into sequence record objects
	for Sixteen in SixTeens:
		SixteenFasta = getDNAFasta(Sixteen)
		SixteenSSubunits.append(SixteenFasta)
	handle2.close()
#-------------------------------------------------------------------------------------------------
# 5: Converts sequence record object as a reverse complement FASTA formated sequnece.
def getReverseComplementFasta(SeqRecord):
	reverseCompSeq = SeqRecord.seq.reverse_complement()
	SeqRecord.seq = reverseCompSeq
	FASTA = SeqRecord.format("fasta") 
	FASTA = fastaClean(FASTA)
	return FASTA
#-------------------------------------------------------------------------------------------------
# 6: Cleans up FASTA formated sequences.
def fastaClean(FASTA):
	FASTAHeader, FASTACode = FASTA.split("\n", 1) # Splits FASTA's into header and genetic code.
	FASTACode = FASTACode.replace("-","").replace("\n","") # Removes alignment markers and converts FASTA file sequence into a single line.
	FASTA = FASTAHeader + "\n" + FASTACode
	return FASTA
#-------------------------------------------------------------------------------------------------
# 7: Creates a more informative header for the 16S gene.
def fastaHeaderSwap(FASTA, subjectAccession):
	FASTAHeader, FASTACode = FASTA.split("\n", 1) # Splits FASTA's into header and genetic code.
	FASTAHeader = ">" + subjectAccession + " 16S rRNA gene"
	FASTA = FASTAHeader + "\n" + FASTACode
	return FASTA
#-------------------------------------------------------------------------------------------------
# 8: Appends genome accession to a file that acts as a list of bad accessions..
def appendBadGenomeList(genome):
	badAccession = path.split(genome)[1].strip(".fna")
	try:
		oufile = open("No16SGenomesHMM.txt", "a")
		oufile.write(badAccession + "\n")
		oufile.close()
	except IOError:
		print "Failed to open " + oufile
		exit(1)
#-------------------------------------------------------------------------------------------------
# 9: Adds SixteenS gene to a FASTA file.
def write16SToFile(SixteenSGene):
	try:
		oufile = open("Found16SGenesHMM.fna", "a")
		oufile.write(SixteenSGene + "\n")
		oufile.close()
	except IOError:
		print "Failed to open " + oufile
		exit(1)
#===========================================================================================================
# Main program code:
# House keeping...
argsCheck(2) # Checks if the number of arguments are correct.

genome = sys.argv[1]
print "Opening " + genome + "..."

subjectAccession = path.split(genome)[1].strip(".fna")

# File extension check
if not genome.endswith(".fna"):
	print "[Warning] " + genome + " may not be a nucleic acid fasta file!"
	
HMMERDBFile = sys.argv[2]
print "Opening " + HMMERDBFile + "..."
print "Searching " + genome + " with " + HMMERDBFile + "..."

SixteenSSubunits = [] 
Found16S = False
try:
	inFile = open(genome, "rU")
	FASTA = inFile.read()
	inFile.close()
	Found16S = runHMMSearch(FASTA, HMMERDBFile) # Pass this FASTA to hmmsearch. runHMMSearch returns true if a 16S was found.
	if Found16S == True: # If we get a result from hmmsearch, check the alignment file.
		print "Found a 16S in the positive strand."
		add16SSequences(SixteenSSubunits)
	else:
		print "No 16S found in the positive strand."
	
	handle = cStringIO.StringIO(FASTA) # Instead reading from the file again we make a virtual file from a string and pass this.
	FASTA = []
	SeqRecords = SeqIO.parse(handle, "fasta")	
	for record in SeqRecords:
		FASTA.append(getReverseComplementFasta(record))
	FASTA = "\n".join(FASTA)
	handle.close()

	Found16S = runHMMSearch(FASTA, HMMERDBFile) # Pass this FASTA to hmmsearch. runHMMSearch returns true if a 16S was found.
	if Found16S == True: # If we get a result from hmmsearch, check the alignment file.
		print "Found a 16S in the negitive strand."
		add16SSequences(SixteenSSubunits)
	else:
		print "No 16S found in the negitive strand."
except IOError:
	print "Failed to open " + inFile
	exit(1)
	
Top16S = ""
if SixteenSSubunits:
	
	Top16SLength = 0
	for s in SixteenSSubunits:
		Current16SSeqLength = len(s.split("\n", 1)[1]) # Splits FASTA into a two element list (Header, Sequence). Gets length of the Seq.
		if Current16SSeqLength >= Top16SLength:
			Top16S = s
			Top16SLength = Current16SSeqLength
	if len(Top16S) < 2000 and len(Top16S) > 1000: # 16S genes are around 1500 B.P. This fitres out partial sequence or really large sequences.
		Top16S = fastaHeaderSwap(Top16S, subjectAccession)
		write16SToFile(Top16S)
		print "Writing best 16S to file."
	else:
		appendBadGenomeList(genome) # If 16S gene is too partial to be used.
		print "Though a partial 16S was found, it was of low quality." 
		print "Writing genome accession to No16SGenomesHMM.txt"  
else:
	appendBadGenomeList(genome)
	print "No 16S found. Writing genome accession to No16SGenomesHMM.txt" 
print "Done!\n"