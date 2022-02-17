#!/usr/bin/env python3
"""
Author : eijy <eijy@localhost>
Date   : 2022-02-15
Purpose: Tretranucleotide frequency
"""

import argparse
import os
from typing import NamedTuple, Dict
from collections import Counter

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

    counts = Counter(args.dna)

    print(counts.get('A',0), 
          counts.get('C',0), 
          counts.get('G',0), 
          counts.get('T',0))



# --------------------------------------------------
def test_count() -> None:
    """ Test count """

    assert Counter('') == {}
    assert Counter('123XYZ') == {'1':1, '2':1, '3':1, 'X':1, 'Y':1, 'Z':1}
    assert Counter('A') == {'A':1}
    assert Counter('C') == {'C':1}
    assert Counter('G') == {'G':1}
    assert Counter('T') == {'T':1}
    assert Counter('ACCGGGTTTT') == {'A':1, 'C':2, 'G':3, 'T':4}



# --------------------------------------------------
if __name__ == '__main__':
    main()
