#!/usr/bin/env bash

filename="$1"
samtools faidx $filename #in fasta format
cut -f1,2 ${filename}.fai > sizes.genome
