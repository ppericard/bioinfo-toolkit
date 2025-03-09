#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Split Paired Fastq

Description : Split a paired fastq into 2 files
  
  # Running examples:
  
  split_paired_fastq.py -i input.fastq
  # output files are input.left_reads.fastq and input.right_reads.fastq
  
  split_paired_fastq.py -i input.fastq \
                        -1 output.left_reads.fastq \
                        -2 output.right_reads.fastq

  split_paired_fastq.py < input.fastq
  split_paired_fastq.py -i input.fastq
  # output files are left_reads.fastq and right_reads.fastq

-----------------------------------------------------------------------

Author: This software is written and maintained by Pierre Pericard
(pierre.pericard@ed.univ-lille1.fr)
Created: 2016-04-13
Last Modified: 2016-04-13
Licence: GNU GPL 3.0

Copyright 2016 Pierre Pericard

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


def read_fastq_file_handle(fastq_file_handle):
    """
    Parse a fastq file and return a generator
    """
    
    # Variables initialization
    line_count = 0
    header = ''
    seq = ''
    qual = ''
    
    # Reading input file
    for line in (l.strip() for l in fastq_file_handle if l.strip()):
        line_count += 1
        if line_count % 4 == 1:
            if header:
                # Yield previous sequence
                yield header, seq, qual
            # Get read complete header (sequence id and description)
            header = line[1:]
        elif line_count % 4 == 2:
            # Get read sequence
            seq = line
        elif line_count % 4 == 0:
            # Get read quality
            qual = line
    
    # Yield the last sequence
    yield header, seq, qual
    
    # Close input file
    fastq_file_handle.close()


if __name__ == '__main__':

    # Argument parser initialization
    parser = argparse.ArgumentParser(description='Split a paired fastq into 2 files')
    
    # -i / --input_fastq
    parser.add_argument('-i', '--input_fastq',
                        action='store',
                        metavar='INFASTQ', 
                        type=argparse.FileType('r'),
                        default='-',
                        help="Input paired fastq file. "
                             "Default to <stdin>")
    
    # -1 / --output_left
    parser.add_argument('-1', '--output_left', 
                        action='store',
                        metavar='OUTLEFT', 
                        type=str, 
                        help="Output left fastq file.")
    
    # -2 / --output_right
    parser.add_argument('-2', '--output_right',
                        action='store',
                        metavar='OUTRIGHT', 
                        type=str, 
                        help="Output right fastq file.")
    
    # --fixed_names # Hidden argument
    parser.add_argument('--fixed_names',
                        action='store_true',
                        help=argparse.SUPPRESS)
    
    # Parse arguments from command line
    args = parser.parse_args()
    
    # Set fixed names if input fastq file is read from STDIN
    if args.input_fastq.name == '<stdin>':
        args.fixed_names = True
    
    # Set default output filenames if none is provided,
    # otherwise output reads only to the provided files
    # (either left, right, or both)
    if not (args.output_left or args.output_left):
        input_fastq_basename = '.'.join(args.input_fastq.name.split('.')[:-1])
        args.output_left = input_fastq_basename + '.left_reads.fastq'
        args.output_right = input_fastq_basename + '.right_reads.fastq'
    
    # Use fixed filenames if set.
    # Is needed by some workflow managers (eg. Galaxy)
    if args.fixed_names:
        args.output_left = 'left_reads.fastq'
        args.output_right = 'right_reads.fastq'
    
    # Open right and left output files
    output_left_fh = None
    output_right_fh = None
    if args.output_left:
        output_left_fh = open(args.output_left, 'w')
    if args.output_right:
        output_right_fh = open(args.output_right, 'w')
    
    # Variables initialization
    seq_count = 0
    
    # Read input fastq file:
    for header, seq, qual in read_fastq_file_handle(args.input_fastq):
        seq_count += 1
        if (seq_count % 2 == 1): # is left read
            if output_left_fh:
                output_left_fh.write('@{0}\n{1}\n+\n{2}\n'.format(header, seq, qual))
        else: # is right read
            if output_right_fh:
                output_right_fh.write('@{0}\n{1}\n+\n{2}\n'.format(header, seq, qual))
    
