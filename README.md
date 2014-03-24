Phylogenetic-Tree-Building
==========================
A tool kit for extracting 16S rRNA genes from genomes for the purpose of building phylogenetic trees.

* **16SBLAST.py** - Uses NCBI's BLASTn to search for 16S genes within a query genome by BLASTing towards a BLAST database which contains a verity of 16S genes.
* **16SHMMER.py** - Uses [HMMER](http://hmmer.janelia.org) and a 16S hmm to search for 16S genes in a target genome. Searchs both the forward and reverse strand of the genome for the best 16S gene.
* **runPrimerSearch16S.sh** - A shell script that uses the command-line tool grep to search for 16S primer binding sites within a genome. 
* **runHMMSearchGyrase.sh** - A shell sript that uses [HMMER](http://hmmer.janelia.org) to search for Gyrase B. genes within a genome.
