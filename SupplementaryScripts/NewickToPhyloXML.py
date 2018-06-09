#!/usr/bin/env python 
# Created by: Lee Bergstrand 
# Description: A simple program that converts Newick phylogenetic tree format to the newer PhyloXML format.
#
# Requirements: - This script requires the Biopython module: http://biopython.org/wiki/Download
#
# Usage: NewickToPhyloXML.py <PhyloTree.nwk>
# Example: NewickToPhyloXML.py PhyloTree.nwk
# ----------------------------------------------------------------------------------------
# ===========================================================================================================

# Imports:

import sys

from Bio import Phylo

# ===========================================================================================================
# Functions:

# 1: Checks if in proper number of arguments are passed gives instructions on proper use.
def argsCheck(numArgs):
    if len(sys.argv) < numArgs or len(sys.argv) > numArgs:
        print("Sequence Downloader")
        print("By Lee Bergstrand\n")
        print("Usage: " + sys.argv[0] + " <PhyloTree.nwk>")
        print("Examples: " + sys.argv[0] + " PhyloTree.nwk\n")
        exit(1)  # Aborts program. (exit(1) indicates that an error occurred)


# ===========================================================================================================
# Main program code:

# House keeping...
argsCheck(2)  # Checks if the number of arguments are correct.

# Stores file one for input checking.
print(">> Opening Newicktree...")
inFile = sys.argv[1]
outFile = inFile.split(".")[0] + ".xml"

# File extension check
if not inFile.endswith(".nwk"):
    print("[Warning] " + inFile + " may not be a Newick file!")

print(">> Converting to PhyloXML...")
# Converts Newick to PhyloXML.
try:
    Phylo.convert(inFile, 'newick', outFile, 'phyloxml')
except IOError:
    print("Failed to open " + inFile + "or" + outFile)
    exit(1)

print(">> Done...")
