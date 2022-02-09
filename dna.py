#!/usr/bin/env python3
"""
Author : eijy <eijy@localhost>
Date   : 2022-01-12
Purpose: Tretranucleotide frequency
"""

import argparse
import os
from typing import NamedTuple, Tuple

class Args(NamedTuple):
    """ Command-line arguments """
    dna: str

# --------------------------------------------------
def get_args() -> Args:
    """ Get command-line arguments """

    parser = argparse.ArgumentParser(
        description='Tretranucleotide frequency',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('dna',
                        metavar='DNA',
                        help='Input DNA sequence')

    args = parser.parse_args()

    if os.path.isfile(args.dna):
        args.dna = open(args.dna).read().rstrip()

    return Args(args.dna)


# --------------------------------------------------
def main() -> None:
    """ Count nucleotide frequency """
    args = get_args()

    count_a, count_c, count_g, count_t = count(args.dna)

    #print('{} {} {} {}'.format(count_a, count_c, count_g, count_t))

    print(f'{count_a} {count_c} {count_g} {count_t}')


# --------------------------------------------------
def count(dna: str) -> Tuple[int, int, int, int]:
    """ Count bases in DNA """

    #count_a, count_c, count_g, count_t = 0, 0, 0, 0
    #for base in dna:
    #    if base == 'A':
    #        count_a += 1
    #    elif base == 'C':
    #        count_c += 1
    #    elif base == 'G':
    #        count_g += 1
    #    elif base == 'T':
    #        count_t += 1
    #return(count_a, count_c, count_g, count_t)

    # approach2
    return (dna.count('A'), dna.count('C'), dna.count('G'), dna.count('T'))





# --------------------------------------------------
if __name__ == '__main__':
    main()
