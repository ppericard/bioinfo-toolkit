#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
"""

import argparse
import pandas




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
    for i in xrange(0, len(seq), linereturn):
        buff.append("{0}\n".format(seq[i:(i + linereturn)]))
    return ''.join(buff).rstrip()


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Filter a fasta file based on sequence length.')
    parser.add_argument('input_fasta_files', metavar='FASTA', 
                        type=argparse.FileType('r'), nargs='+',
                        help='One or more fasta files')
    parser.add_argument('-o', '--output_histo', metavar='PDF', 
                        type=argparse.FileType('wb'), default='histo.pdf',
                        help='ouput pdf histogram file')
    parser.add_argument('--x_log', action='store_true',
                        help='X axis log scale')
    parser.add_argument('--y_log', action='store_true',
                        help='Y axis log scale')
    parser.add_argument('-b', '--y_bin', metavar='INT',
                        type=int, default=10,
                        help='Y axis bin size. '
                             'Default is %(default)s')
    args = parser.parse_args()
    
    # Code is duplicated here to prevent to have to test args.max_length many times
    if args.max_length:
        for header, sequence in read_fasta_file_handle(args.input_fasta):
            if len(sequence) >= args.min_length and len(sequence) <= args.max_length:
                args.output_fasta.write(">{0}\n{1}\n".format(header, format_seq(sequence)))
    else:
        for header, sequence in read_fasta_file_handle(args.input_fasta):
            if len(sequence) >= args.min_length:
                args.output_fasta.write(">{0}\n{1}\n".format(header, format_seq(sequence)))
