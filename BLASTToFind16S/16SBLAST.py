#!/usr/bin/env python

# ----------------------------------------------------------------------------------------
# Created by: Lee Bergstrand
# Description: A simple python program that use BLASTn to search for 16S genes within a geneome.
#              Extracts the aligned query sequence from the BLAST results (This should be the 16S
#              gene in the query genome).
#
# Requirements: - This program requires the Biopython module: http://biopython.org/wiki/Download
#               - This script requires BLAST+ 2.2.9 or later.
#               - MakeNABlastDB must be used to create BLASTn databases for both query and subject proteomes.
#               - BLAST databases require that the FASTA file they were made from remain in the same directory.
#  
# Usage: 16SBLAST.py <QueryGenome.faa> <16SDataBase.faa>
# Example: 16SBLAST.py QueryGenome.faa RDPActinoBacter16S.faa
# ----------------------------------------------------------------------------------------
# ===========================================================================================================

# Imports & Setup:
import csv
import subprocess
import sys
from multiprocessing import cpu_count
from os import path

processors = cpu_count()  # Gets number of processor cores for BLAST.


# ===========================================================================================================
# Functions:

# 1: Checks if in proper number of arguments are passed gives instructions on proper use.
def argsCheck():
    if len(sys.argv) < 3:
        print("16S Gene Finder")
        print("By Lee Bergstrand\n")
        print("Please refer to source code for documentation\n")
        print("Usage: " + sys.argv[0] + "<QueryGenome.faa> <16SDataBase.faa>\n")
        print("Examples:" + sys.argv[0] + "QueryGenome.faa RDPActinoBacter16S.faa")
        exit(1)  # Aborts program. (exit(1) indicates that an error occurred)


# -------------------------------------------------------------------------------------------------
# 2: Runs BLASTn with settings specific for extracting subject sequences.
def runBLASTFor16S(query, BLASTDBFile):
    # Runs BLASTn and saves the output to a string.
    # Blastn is set to output a csv which can be parsed by Pythons CSV module.
    BLASTOut = subprocess.check_output(
        ["blastn", "-db", BLASTDBFile, "-query", query, "-num_threads", str(processors), "-outfmt",
         "10 qseqid sseqid length qseq evalue bitscore"])
    return BLASTOut


# -------------------------------------------------------------------------------------------------
# 3: Appends genome accession to a file that acts as a list of bad accessions.
def appendBadGenomeList(genome):
    global outfile
    badAccession = path.split(genome)[1].strip(".fna")
    try:
        outfile = open("No16SGenomesBLAST.txt", "a")
        outfile.write(badAccession + "\n")
        outfile.close()
    except IOError:
        print("Failed to open {0}".format(outfile))
        exit(1)


# -------------------------------------------------------------------------------------------------
# 4: Cleans up FASTA formatted sequences.
def fastaClean(FASTA):
    FASTAHeader, FASTACode = FASTA.split("\n", 1)  # Splits FASTA's into header and genetic code.
    # Removes alignment markers and converts FASTA file sequence into a single line.
    FASTACode = FASTACode.replace("-", "").replace("\n", "")
    FASTA = FASTAHeader + "\n" + FASTACode
    return FASTA


# ===========================================================================================================
# Main program code:
# House keeping...
argsCheck()  # Checks if the number of arguments are correct.

queryFile = sys.argv[1]
print("Opening " + queryFile + "...")

subjectAccession = path.split(queryFile)[1].strip(".fna")

# File extension check
if not queryFile.endswith(".fna"):
    print("[Warning] " + queryFile + " may not be a nucleic acid fasta file!")

BLASTDBFile = sys.argv[2]
print("Opening " + BLASTDBFile + "...")
print("Blasting " + queryFile + " against " + BLASTDBFile + "...")
BLASTResults = runBLASTFor16S(queryFile, BLASTDBFile)

if not BLASTResults:  # If there are no BLAST results aborts the program.
    print("\nUnfortunately there are no 16S BLAST results for " + queryFile + ". Try using another BLAST DB. or")
    print("you may also want to try another method to find 16S other than BLAST (eg. HMMs).\n")
    print("Writing genome accession to No16SGenomesBLAST.txt")
    appendBadGenomeList(queryFile)
    exit(1)  # Aborts program. (exit(1) indicates that an error occurred)

BLASTCSVOut = BLASTResults.splitlines(True)  # Converts raw BLAST csv output into list of csv rows.
BLASTreader = csv.reader(BLASTCSVOut)  # Reads BLAST csv rows as a csv.

Found16S = False
Top16S = ""
Top16SLength = 0
for row in BLASTreader:
    Current16S = row[3]
    Current16SLength = len(Current16S)
    # 16S genes are around 1500 B.P. Below filters out partial sequence or really large sequences.
    if 2000 > Current16SLength > 1000:
        if Current16SLength > Top16SLength:
            Found16S = True
            Top16S = Current16S
            Top16SLength = Current16SLength

if not Found16S:  # If there are no BLAST that are around the size of 16S rRNA abort the program.
    print("\nUnfortunately there are no 16S BLAST results for " + queryFile + ". Try using another BLAST DB or")
    print("you may also want to try another method to find 16S other than BLAST (eg. HMMs).\n")
    print("Writing genome accession to No16SGenomesBLAST.txt")
    appendBadGenomeList(queryFile)
    exit(1)  # Aborts program. (exit(1) indicates that an error occurred)

print("Extracting 16S BLAST Results!")
FASTA = ">" + subjectAccession + "\n" + Top16S
FASTA = fastaClean(FASTA)

print("Writing results to file.")
try:
    fileWriter = open("Found16SGenesBLAST.fna", "a")
    fileWriter.write(FASTA + "\n")
    fileWriter.close()
except IOError:
    print("Error writing " + "Found16SGenesBLAST.fna" + " to file.")
print("Done.\n")
