#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
GenerSampleFasta

Description: Generate a sample given a fasta file and a sample size
               or a draw probability

  generSampleFasta.py -i input.fa -o sample.fa -n 4200
  generSampleFasta.py -i input.fa -o sample.fa -p 0.01

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

    random.seed(os.urandom(128))

    # Arguments
    parser = argparse.ArgumentParser(description='',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input_fasta', metavar='input',
                        type=argparse.FileType('r'), required=True,
                        help='input fasta file')
    parser.add_argument('-o', '--output_sample', metavar='sample',
                        type=argparse.FileType('w'), required=True,
                        help='ouput sample fasta file')
    parser.add_argument('--output_non_sample', metavar='nonsample',
                        type=argparse.FileType('w'), required=False,
                        help='ouput non sample fasta file')
    parser.add_argument('-p', '--proba', metavar='proba',
                        type=float, default=0.01,
                        help='Probability to keep a sequence')
    parser.add_argument('-n', '--seq_nb', metavar='number',
                        type=int, help='Nb of sequences to keep')
    args = parser.parse_args()

    if args.seq_nb:
        input_seq_nb = int(subprocess.check_output("grep -c '>' {0}".format(args.input_fasta.name), shell=True).strip())
        bool_list = [False for i in range(input_seq_nb)]
        for random_index in random.sample([i for i in range(input_seq_nb)], args.seq_nb):
            bool_list[random_index] = True
        for header, sequence in read_fasta_file_handle(args.input_fasta):
            if bool_list.pop():
                args.output_sample.write(">{0}\n{1}\n".format(header, format_seq(sequence)))
            else:
                if args.output_non_sample:
                    args.output_non_sample.write(">{0}\n{1}\n".format(header, format_seq(sequence)))
    else:
        for header, sequence in read_fasta_file_handle(args.input_fasta):
            if random.random() <= args.proba:
                args.output_sample.write(">{0}\n{1}\n".format(header, format_seq(sequence)))
            else:
                if args.output_non_sample:
                    args.output_non_sample.write(">{0}\n{1}\n".format(header, format_seq(sequence)))
