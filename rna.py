#!/usr/bin/env python3
"""
Author : eijy <eijy@localhost>
Date   : 2022-02-22
Purpose: Transcribe DNA to RNA
"""

import argparse
from typing import NamedTuple, TextIO


class Args(NamedTuple):
    """ Command-line arguments """
    positional: str
    string_arg: str
    int_arg: int
    file: TextIO
    on: bool


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
                        types=argparse.FileType('rt'))


    args = parser.parse_args()

    return Args(args.positional, args.arg, args.int, args.file, args.on)


# --------------------------------------------------
def main() -> None:
    """ Make a jazz noise here """

    args = get_args()
    str_arg = args.string_arg
    int_arg = args.int_arg
    file_arg = args.file
    flag_arg = args.on
    pos_arg = args.positional

    print(f'str_arg = "{str_arg}"')
    print(f'int_arg = "{int_arg}"')
    print('file_arg = "{}"'.format(file_arg.name if file_arg else ''))
    print(f'flag_arg = "{flag_arg}"')
    print(f'positional = "{pos_arg}"')


# --------------------------------------------------
if __name__ == '__main__':
    main()
