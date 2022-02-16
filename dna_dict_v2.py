#!/usr/bin/env python3
"""
Author : eijy <eijy@localhost>
Date   : 2022-02-15
Purpose: Tretranucleotide frequency
"""

import argparse
import os
from typing import NamedTuple, Tuple
from collections import defaultdict

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

    counts = count(args.dna)

    print(counts.get('A',0), 
          counts.get('C',0), 
          counts.get('G',0), 
          counts.get('T',0))



# --------------------------------------------------
def count(dna: str) -> Dict[str, int]:
    """ Count bases in DNA """

    counts: Dict[str, int] = defaultdict(int)

    for base in dna:
        counts[base] += 1

    return counts



# --------------------------------------------------
def test_count() -> None:
    """ Test count """

    assert count('') == {}
    assert count('123XYZ') == {'1':1, '2':1, '3':1, 'X':1, 'Y':1, 'Z':1}
    assert count('A') == {'A':1}
    assert count('C') == {'C':1}
    assert count('G') == {'G':1}
    assert count('T') == {'T':1}
    assert count('ACCGGGTTTT') == {'A':1, 'C':2, 'G':3, 'T':4}



# --------------------------------------------------
if __name__ == '__main__':
    main()
