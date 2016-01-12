#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
GenerSampleFastq

Description: Generate a sample given a fastq file and a sample size
               or a draw probability

  generSampleFastq.py -i input.fq -o sample.fq -n 4200
  generSampleFastq.py -i input.fq -o sample.fq -p 0.01

-----------------------------------------------------------------------

Author: This software is written and maintained by Pierre Pericard
(pierre.pericard@ed.univ-lille1.fr)
Created: 2014
Last Modified: 2016-01-11
Licence: GNU GPL 3.0

Copyright 2014-2016 Pierre Pericard

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
import random
import argparse
import subprocess

def read_fastq_file_handle(fastq_file_handle):
    """
    Parse a fastq file and return a generator
    """
    # Variables initialization
    header = ''
    sequence = ''
    quality = ''
    line_count = 0
    seq_count = 0
    # Reading input file
    for line in fastq_file_handle:
        line = line.strip()
        if line:
            line_count += 1
            # Yield the last read header, sequence and quality
            if line_count % 4 == 1:
                if seq_count:
                    yield (header, sequence, quality)
                # Get header
                header = line[1:]
                seq_count += 1
            elif line_count % 4 == 2:
                sequence = line
            elif line_count % 4 == 0:
                quality = line
    # Yield the input file last read
    yield (header, sequence, quality)
    # Close input file
    fastq_file_handle.close()

if __name__ == '__main__':
    
    random.seed(os.urandom(128))
    
    # Arguments
    parser = argparse.ArgumentParser(description='', 
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input_fastq', metavar='input', 
                        type=argparse.FileType('r'), required=True,
                        help='input fastq file')
    parser.add_argument('-o', '--output_sample', metavar='sample',
                        type=argparse.FileType('w', 0), required=True,
                        help='ouput sample fastq file')
    parser.add_argument('--output_non_sample', metavar='nonsample',
                        type=argparse.FileType('w', 0), required=False,
                        help='ouput non sample fastq file')
    parser.add_argument('-p', '--proba', metavar='proba',
                        type=float, default=0.01,
                        help='Probability to keep a sequence')
    parser.add_argument('-n', '--seq_nb', metavar='number',
                        type=int, help='Nb of sequences to keep')
    args = parser.parse_args()
    
    line_count = 0
    seq_count = 0
    
    if args.seq_nb:
        input_seq_nb = int(subprocess.check_output("wc -l {0}".format(args.input_fastq.name), shell=True).split()[0].strip()) / 4
        bool_list = [False for i in range(input_seq_nb)]
        for random_index in random.sample([i for i in range(input_seq_nb)], args.seq_nb):
            bool_list[random_index] = True
        for header, sequence, quality in read_fastq_file_handle(args.input_fastq):
            read_string = "@{0}\n{1}\n+\n{2}\n".format(header, sequence, quality)
            if bool_list.pop():
                args.output_sample.write(read_string)
            else:
                if args.output_non_sample:
                    args.output_non_sample.write(read_string)
    else:
        for header, sequence, quality in read_fastq_file_handle(args.input_fastq):
            read_string = "@{0}\n{1}\n+\n{2}\n".format(header, sequence, quality)
            if random.random() <= args.proba:
                args.output_sample.write(read_string)
            else:
                if args.output_non_sample:
                    args.output_non_sample.write(read_string)
