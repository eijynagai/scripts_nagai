#!/usr/bin/env python3
"""
Author : eijy <eijy@localhost>
Date   : 2022-02-28
Purpose: Transcribe DNA to RNA
"""

import argparse
from typing import NamedTuple, TextIO


class Args(NamedTuple):
    """ Command-line arguments """
    files: List[TextIO]
    out_dir: str

# --------------------------------------------------
def get_args() -> Args:
    """ Get command-line arguments """

    parser = argparse.ArgumentParser(
        description='Transcribe DNA to RNA',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-o',
                        '--out_dir',
                        help='Outpud directory',
                        metavar='DIR',
                        type=str,
                        default='out')

    parser.add_argument('file',
                        help='Input DNA file(s)',
                        metavar='FILE',
                        nargs='+',
                        type=argparse.FileType('rt'))


    args = parser.parse_args()

    return Args(files=args.file, out_dir=args.out_dir)


# --------------------------------------------------
def main() -> None:
    """ Make a jazz noise here """

    args = get_args()
    
    if not os.path.isdir(args.out_dir):
        os.makedirs(args.out_dir)

    num_files, num_seqs = 0, 0
    for fh in args.files:
    # open an output file in the output directory
    
    # for each line/sequence from the input file:
        # write the transcribed sequence to the output file
        # update the number of files processed
    print('Done!')

# --------------------------------------------------
if __name__ == '__main__':
    main()
