#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Fasta Name Filter

Description: Filter a fasta file based on a string to find in the
               sequences headers, or given a file with a list of id

  fasta_name_filter.py -i input.fa -o output.fa -p "pattern"
  fasta_name_filter.py -i input.fa -o output.fa -f patterns.ids

-----------------------------------------------------------------------

Author: This software is written and maintained by Pierre Pericard
(pierre.pericard@gmail.com)
Created: 2014
Last Modified: 2017-02-12
Licence: GNU GPL 3.0

Copyright 2014-2017 Pierre Pericard

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

import sys
import os
import argparse
import string
import re

def read_fasta_file_handle(fasta_file_handle):
    """
    Parse a fasta file and return a generator
    """
    # Variables initialization
    header = ''
    seqlines = list()
    sequence_nb = 0
    # Reading input file
    for line in fasta_file_handle:
        if line[0] == '>':
            # Yield the last read header and sequence
            if sequence_nb:
                yield (header, ''.join(seqlines))
                del seqlines[:]
            # Get header
            header = line[1:].rstrip()
            sequence_nb += 1
        else:
            # Concatenate sequence
            seqlines.append(line.strip())
    # Yield the input file last sequence
    yield (header, ''.join(seqlines))
    # Close input file
    fasta_file_handle.close()

def format_seq(seq, linereturn=80):
    """
    Format an input sequence
    """
    buff = list()
    for i in range(0, len(seq), linereturn):
        buff.append("{0}\n".format(seq[i:(i + linereturn)]))
    return ''.join(buff).rstrip()


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = 'Filter a fasta file based on pattern in the sequence name.')
    # -i / --input_fasta
    parser.add_argument('-i', '--input_fasta',
                        action  = 'store',
                        metavar = 'INFILE',
                        type    = argparse.FileType('r'),
                        default = '-',
                        help    = 'Input fasta file')
    # -o / --output_accepted
    parser.add_argument('-o', '--output_accepted',
                        action   = 'store',
                        metavar  = 'OUTFILE',
                        type     = argparse.FileType('w'),
                        required = False,
                        help     = 'Output accepted fasta file')
    # --output_rejected
    parser.add_argument('--output_rejected',
                        action   = 'store',
                        metavar  = 'OUTFILE',
                        type     = argparse.FileType('w'),
                        required = False,
                        help     = 'Output rejected fasta file')
    # -p / --pattern
    parser.add_argument('-p', '--pattern',
                        action  = 'store',
                        metavar = 'STR',
                        type    = str,
                        help    = 'Pattern to find')
    # -f / --patterns_file
    parser.add_argument('-f', '--patterns_file',
                        action  = 'store',
                        metavar = 'INFILE',
                        type    = argparse.FileType('r'),
                        help    = 'File with patterns to find')
    # --ids
    parser.add_argument('--ids',
                        action = 'store_true',
                        help   = 'Pattern(s) is/are sequence ids. '
                                 'Default is to search the pattern in the complete header')
    #
    args = parser.parse_args()

    # Check that at least a pattern or a file is provided
    if not args.pattern and not args.patterns_file:
        parser.print_help()
        raise Exception('A pattern or a patterns file has to be supplied')

    # If neither accepted nor rejected output filepath is provided
    # revert to outputing accepted sequences to STDOUT
    if not (args.output_accepted or args.output_rejected):
        args.output_accepted = sys.stdout

    # Store patterns
    patterns_list = list()
    if args.patterns_file:
        # read patterns and store them
        for pattern in (l.rstrip() for l in args.patterns_file if l.strip()):
            patterns_list.append(pattern)
    else:
        patterns_list.append(args.pattern)


    if args.ids:
        # We search for exact ids matches
        # Convert the patterns list to a frozenset for fast search
        patterns_set = frozenset(patterns_list)
        for header, sequence in read_fasta_file_handle(args.input_fasta):
            seq_id = header.split()[0]
            if seq_id in patterns_set:
                if args.output_accepted:
                    args.output_accepted.write(">{0}\n{1}\n".format(header, format_seq(sequence)))
            else:
                if args.output_rejected:
                    args.output_rejected.write(">{0}\n{1}\n".format(header, format_seq(sequence)))
    else:
        # We search for matches in the header
        re_patterns = re.compile('|'.join(patterns_list), flags=re.IGNORECASE)
        for header, sequence in read_fasta_file_handle(args.input_fasta):
            if re_patterns.search(header):
                if args.output_accepted:
                    args.output_accepted.write(">{0}\n{1}\n".format(header, format_seq(sequence)))
            else:
                if args.output_rejected:
                    args.output_rejected.write(">{0}\n{1}\n".format(header, format_seq(sequence)))
