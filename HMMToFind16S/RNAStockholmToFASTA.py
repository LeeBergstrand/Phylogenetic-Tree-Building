#!/usr/bin/env python 
# Created by: Lee Bergstrand 
# Descript: A simple program that takes a RNA nucleotide stockholm file and returns a DNA FASTA file.
#
# Requirements: - This script requires the Biopython module: http://biopython.org/wiki/Download
#     
# Usage: getGenbankSeqs.py <sequences.sto> 
# Example: getGenbankSeqs.py mySeqs.sto
# ----------------------------------------------------------------------------------------
# ===========================================================================================================

# Imports:
import sys

from Bio import SeqIO

# ===========================================================================================================
# Functions:

# 1: Checks if in proper number of arguments are passed gives instructions on proper use.
def argsCheck(numArgs):
    if len(sys.argv) < numArgs or len(sys.argv) > numArgs:
        print("A simple program that takes a RNA nucleotide stockholm file and returns a DNA FASTA file.")
        print("By Lee Bergstrand\n")
        print("Usage: " + sys.argv[0] + " <sequences.txt> [email@mail.com]")
        print("Examples: " + sys.argv[0] + " mySeqs.txt JBro@YOLO.com\n")
        exit(1)  # Aborts program. (exit(1) indicates that an error occurred)


# 2: Converts the sequence record object from an RNA Stockholm file to a DNA FASTA file.
def getDNAFasta(SeqRecord):
    SeqRecord.letter_annotations = {}  # Removes per letter annitations. Biopython throws an error if you try to
    # reverse complement a sequence with per letter annotations (since afterward
    # these annotations would no longer be valid). We strip these per letter
    # annotions since they are not a part of the FASTA format anyways.
    DNA = SeqRecord.seq.back_transcribe()
    SeqRecord.seq = DNA
    return SeqRecord.format("fasta")


# ===========================================================================================================
# Main program code:

# House keeping...
argsCheck(2)  # Checks if the number of arguments are correct.

# Stores file one for input checking.
print(">> Opening Stockholm file...")
inFile = sys.argv[1]
outFile = inFile + ".fna"

# File extension check
if not inFile.endswith(".sto"):
    print("[Warning] " + inFile + " may not be a Stockholm file!")

try:
    writer = open(outFile, "w")
    handle = open(inFile, "rU")
    SeqRecords = SeqIO.parse(handle, "stockholm")
    print(">> Converting to FASTA...")
    for record in SeqRecords:
        FASTA = getDNAFasta(record)
        writer.write(FASTA)
    handle.close()
    writer.close()
except IOError:
    print("Failed to open " + inFile + " or " + outFile)
    exit(1)

print(">> Done.")
