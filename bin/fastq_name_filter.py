#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Fastq Name Filter

Description: Filter a fastq file based on a string to find in the
               sequences headers, or given a file with a list of id

  fastq_name_filter.py -i input.fq -o output.fq -p "pattern"
  fastq_name_filter.py -i input.fq -o output.fq -f patterns.ids

-----------------------------------------------------------------------

Author: This software is written and maintained by Pierre Pericard
(pierre.pericard@gmail.com)
Created: 2014
Last Modified: 2017-04-04
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

def read_fastq_file_handle(fastq_file_handle):
    """
    Parse a fastq file and return a generator
    """
    # Variables initialization
    count = 0
    header = ''
    seq = ''
    qual = ''
    # Reading input file
    for line in (l.strip() for l in fastq_file_handle if l.strip()):
        count += 1
        if count % 4 == 1:
            if header:
                yield header, seq, qual
            header = line[1:]
        elif count % 4 == 2:
            seq = line
        elif count % 4 == 0:
            qual = line
    yield header, seq, qual
    # Close input file
    fastq_file_handle.close()


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = 'Filter a fastq file based on pattern in the sequence name.')
    # -i / --input_fastq
    parser.add_argument('-i', '--input_fastq',
                        action  = 'store',
                        metavar = 'INFILE',
                        type    = argparse.FileType('r'),
                        default = '-',
                        help    = 'Input fastq file')
    # -o / --output_accepted
    parser.add_argument('-o', '--output_accepted',
                        action   = 'store',
                        metavar  = 'OUTFILE',
                        type     = argparse.FileType('w'),
                        required = False,
                        help     = 'Output accepted fastq file')
    # --output_rejected
    parser.add_argument('--output_rejected',
                        action   = 'store',
                        metavar  = 'OUTFILE',
                        type     = argparse.FileType('w'),
                        required = False,
                        help     = 'Output rejected fastq file')
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
        for header, seq, qual in read_fastq_file_handle(args.input_fastq):
            seq_id = header.split()[0]
            if seq_id in patterns_set:
                if args.output_accepted:
                    args.output_accepted.write('@{0}\n{1}\n+\n{2}\n'.format(header, seq, qual))
            else:
                if args.output_rejected:
                    args.output_rejected.write('@{0}\n{1}\n+\n{2}\n'.format(header, seq, qual))
    else:
        # We search for matches in the header
        re_patterns = re.compile('|'.join(patterns_list), flags=re.IGNORECASE)
        for header, seq, qual in read_fastq_file_handle(args.input_fastq):
            if re_patterns.search(header):
                if args.output_accepted:
                    args.output_accepted.write('@{0}\n{1}\n+\n{2}\n'.format(header, seq, qual))
            else:
                if args.output_rejected:
                    args.output_rejected.write('@{0}\n{1}\n+\n{2}\n'.format(header, seq, qual))
