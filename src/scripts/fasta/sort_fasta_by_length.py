#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
sort_fasta_by_length

Description: Sort a FastA file by sequence length (increasing by default)

  sort_fasta_by_length.py --input input.fa --output output.fa

-----------------------------------------------------------------------

Author: This software is written and maintained by Pierre Pericard
(pierre.pericard@univ-lille.fr)
Created: 2015
Last Modified: 2023-04-21
Licence: GNU GPL 3.0

Copyright 2015-2023 Pierre Pericard

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import argparse
import gzip
import sys
from io import TextIOWrapper
from typing import TextIO, Tuple, Iterator
import textwrap


def is_gzip_file(file: TextIO) -> bool:
    file.seek(0)
    magic = file.read(2)
    file.seek(0)
    return magic == b"\x1f\x8b"


def open_input_file(file: TextIO) -> TextIO:
    return (
        TextIOWrapper(gzip.open(file, "rt"))
        if is_gzip_file(file)
        else TextIOWrapper(file)
    )


def parse_fasta_sequences(file: TextIO) -> Iterator[Tuple[str, str]]:
    header, sequence = None, []

    for line in file:
        line = line.rstrip()
        if line.startswith(">"):
            if header:
                yield header, "".join(sequence)
            header, sequence = line, []
        else:
            sequence.append(line)

    if header:
        yield header, "".join(sequence)


def sort_fasta_sequences(
    input_file: TextIO, output_file: TextIO, *, reverse: bool = False
):
    with open_input_file(input_file) as file:
        entries = parse_fasta_sequences(file)

        if not entries:
            raise ValueError("No fasta entries found in the input file.")

        sorted_sequences = sorted(entries, key=lambda x: len(x[1]), reverse=reverse)

    for header, sequence in sorted_sequences:
        print(header, file=output_file)
        print("\n".join(textwrap.wrap(sequence, 80)), file=output_file)


def main():
    parser = argparse.ArgumentParser(
        description="Sort a fasta file by sequence length."
    )
    parser.add_argument(
        "--input",
        type=argparse.FileType("rb"),
        default=sys.stdin.buffer,
        help="Input fasta file, plain text or gzipped. Default: stdin.",
    )
    parser.add_argument(
        "--output",
        type=argparse.FileType("w"),
        default=sys.stdout,
        help="Output sorted fasta file. Default: stdout.",
    )
    parser.add_argument(
        "--reverse",
        action="store_true",
        help="Sort in descending order if set. Default: ascending.",
    )

    args = parser.parse_args()
    sort_fasta_sequences(args.input, args.output, reverse=args.reverse)


if __name__ == "__main__":
    main()

	
