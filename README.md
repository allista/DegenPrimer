DegenPrimer is a tool to check the quality of degenerate primers. 
It calculates melting temperatures, finds stable dimers, hairpins 
and crossdimers, uses BLAST search and full-fledged PCR simulation 
cycle by cycle to provide accurate predictions of possible PCR 
products and their quantities.

All calculations of thermodynamic parameters are made using 
Nearest Neighbor model and the latest thermodynamic tables 
with corrections for particular PCR conditions. PCR simulation 
supports multiplex PCR with arbitrary number of primers.

The program has simple command line interface as well as a 
[separate GUI](https://github.com/allista/DegenPrimerGUI).
PBS is also supported with [PBS-utils](https://github.com/allista/PBS-utils).

Dependecies: 

         python2 >= 2.7
         SciPy >= 0.16
         Biopython >= 1.66
         DendroPy >= 4.0
         BioUtils >= 1.4
