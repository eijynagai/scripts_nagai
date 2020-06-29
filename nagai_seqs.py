#!/usr/bin/env python
# # -*- coding: utf-8 -*-

"""Description
author:     Luis AE Nagai
date:       2020-04-17
Aim:        Sorts of useful scripts for DNA sequences
Input:      DNA sequences
Output:     various

Last modification: 20200417: initial commit
"""

# counting nucleotides from a DNA sequence
def atcg(s):
    return s.count('A'), s.count('T'), s.count('C'), s.count('G')


# transcribing DNA to RNA
def dna2rna(t):
    # replace T to U
    return t.replace('T','U')


# reverse complementary strand
def revcomp(seq):
    complement = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    return "".join(complement[base] for base in seq)[::-1]


# computing GC content of a DNA FASTA file
def gc_content(seq):
    total=len(seq)
    return (seq.count('G')+seq.count('C'))/total*100


def main():
    # runing tests on the functions
    seq='CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGACTGGGAACCTGCGGGCAGTAGGTGGAAT'
    print(gc_content(seq))


    #input fasta file

    #break into ID and sequence

    #calculate gc content

    #add in bigger gc if bigger_gc == ''

    #elif new_gc > bigger_gc, replace it

    #print bigger_gc

if __name__ == '__main__':
    main()
