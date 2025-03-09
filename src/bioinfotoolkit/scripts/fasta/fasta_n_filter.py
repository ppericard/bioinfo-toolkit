#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
fasta_length_filter

Description: Filter a fasta file based on sequence length

  fasta_length_filter.py -i input.fa -o output.fa -m 300
  fasta_length_filter.py -i input.fa -o output.fa -M 1000

-----------------------------------------------------------------------

Author: This software is written and maintained by Pierre Pericard
(pierre.pericard@ed.univ-lille1.fr)
Created: 2019-01-15
Last Modified: 2019-01-15
Licence: GNU GPL 3.0

Copyright 2019 Pierre Pericard

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

import os
import sys
import argparse
import logging
import cProfile
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
    
    # Initiate argument parser
    parser = argparse.ArgumentParser(description='Filter a fasta file based on N content.')
    # -i / --input_fasta
    parser.add_argument('-i', '--input_fasta',
                        action = 'store',
                        metavar = 'FASTA', 
                        type = argparse.FileType('r'),
                        default = '-',
                        help = 'Input fasta file')
    # -o / --output_fasta
    parser.add_argument('-o', '--output_fasta',
                        action = 'store',
                        metavar = 'OUTFASTA', 
                        type=argparse.FileType('w'),
                        default='-',
                        help = 'Ouput fasta file')
    # -r / --max_n_rate
    parser.add_argument('-r', '--max_n_rate', 
                        action = 'store',
                        metavar = 'PCT',
                        type = float,
                        default = 0.1,
                        help='Maximum rate of unknown nucleotides. '
                             'Default is %(default)s')
    # -s / --max_stretch_length
    parser.add_argument('-s', '--max_stretch_length',
                        action = 'store',
                        metavar = 'INT',
                        type = int,
                        default = 5,
                        help = 'Maximum length of allowed streches of unknown nucleotides. '
                               'Default is %(default)s Ns')
    # -a / --replace_n_by_a
    parser.add_argument('-a', '--replace_n_by_a',
                        action = 'store_true',
                        help = 'Replace unknown nucleotides by adenine (N->A)')
    # -t / --dont_trim
    parser.add_argument('-t', '--dont_trim',
                        action = 'store_true',
                        help = 'Do not trim leading and trailing Ns. '
                               'Default is to trim')
    # --debug
    parser.add_argument('--debug',
                        action = 'store_true',
                        help = argparse.SUPPRESS)
    #
    args = parser.parse_args()
    
    # Set logger level
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    
    # Compile Ns homopolymers pattern
    n_regex = re.compile('N+')
    
    #
    for header, sequence in read_fasta_file_handle(args.input_fasta):
        # Trim leading and trailing Ns
        if not args.dont_trim:
            sequence = sequence.strip('nN')
        #
        if sequence.upper().count('N'):
            # Reject sequences with too high a rate of unknown nucleotides
            if sequence.upper().count('N') / len(sequence) > args.max_n_rate:
                logging.debug('REJECT: Above max_n_rate\n>{0}\n{1}\n'.format(header, format_seq(sequence)))
                continue
            #
            len_longest_stretch_n = len(max(n_regex.findall(sequence.upper()), key=len))
            if len_longest_stretch_n > args.max_stretch_length:
                logging.debug('REJECT: Longer than max_stretch_length\n>{0}\n{1}\n'.format(header, format_seq(sequence)))
                continue
            #
            if args.replace_n_by_a:
                sequence = sequence.replace("n", "a").replace("N", "A")
                logging.debug('REPLACE: N->A (>{0})\n'.format(header))
        #
        args.output_fasta.write(">{0}\n{1}\n".format(header, format_seq(sequence)))
