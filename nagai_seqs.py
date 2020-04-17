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



def main():
    # runing tests on the functions
    seq='AAAACCCGGT'
    print(revcomp(seq))



if __name__ == '__main__':
    main()
