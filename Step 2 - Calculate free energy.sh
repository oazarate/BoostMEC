#!/bin/bash

cd data/fasta

RNAfold --noPS -i grna.fa --outfile=grna-free-energy.fa
RNAfold --noPS -i grna-scaffold.fa --outfile=grna-scaffold-free-energy.fa
