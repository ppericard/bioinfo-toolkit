#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
FastQ to FastA

Description: Convert a FastQ file to a FastA file
  
Running examples:
  
  fastq_to_fasta.py -i input.fastq -o output.fasta
  fastq_to_fasta.py < input.fastq > output.fasta

-----------------------------------------------------------------------

Author: This software is written and maintained by Pierre Pericard
(pierre.pericard@ed.univ-lille1.fr)
Created: 2016-04-13
Modified: 2023-06-09
Licence: GNU GPL 3.0

Copyright 2016-2023 Pierre Pericard

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
import sys
import os
from pathlib import Path
from typing import Generator, Tuple, TextIO, Optional, List
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


def read_fastq_file_handle(fastq_file_handle: TextIO) -> Generator[Tuple[str, str, str], None, None]:
    """
    Parse a fastq file and yield sequences as (header, sequence, quality) tuples.
    
    Args:
        fastq_file_handle: An open file handle for the FASTQ file
        
    Yields:
        Tuple containing (header, sequence, quality)
    """
    # Variables initialization
    line_count = 0
    header = ''
    seq = ''
    qual = ''
    
    try:
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
        if header:
            yield header, seq, qual
            
    except Exception as e:
        logger.error(f"Error reading FASTQ file: {e}")
        raise
    finally:
        # Don't close the file if it's stdin
        if fastq_file_handle is not sys.stdin:
            fastq_file_handle.close()


def format_seq(seq: str, line_length: int = 80) -> str:
    """
    Format a sequence with line breaks at specified intervals.
    
    Args:
        seq: The sequence to format
        line_length: Maximum length of each line (default: 80)
        
    Returns:
        Formatted sequence string with line breaks
    """
    return '\n'.join(seq[i:i + line_length] for i in range(0, len(seq), line_length))


def fastq_to_fasta(input_file: TextIO, output_file: TextIO, line_length: int = 80) -> int:
    """
    Convert FASTQ to FASTA format.
    
    Args:
        input_file: Open file handle for input FASTQ
        output_file: Open file handle for output FASTA
        line_length: Length of sequence lines in output (default: 80)
        
    Returns:
        Number of sequences processed
    """
    seq_count = 0
    
    try:
        # Process each sequence from the FASTQ file
        for header, seq, _ in read_fastq_file_handle(input_file):
            output_file.write(f">{header}\n")
            output_file.write(f"{format_seq(seq, line_length)}\n")
            seq_count += 1
            
        return seq_count
    except Exception as e:
        logger.error(f"Error during conversion: {e}")
        return -1


def parse_args() -> argparse.Namespace:
    """Parse and return command line arguments"""
    parser = argparse.ArgumentParser(
        description='Convert a FASTQ file to a FASTA file',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('-i', '--input_fastq',
                        action='store',
                        metavar='INFASTQ', 
                        type=argparse.FileType('r'),
                        default=sys.stdin,
                        help="Input FASTQ file. Default is stdin.")
    
    parser.add_argument('-o', '--output_fasta',
                        action='store',
                        metavar='OUTFASTA', 
                        type=argparse.FileType('w'),
                        default=sys.stdout,
                        help="Output FASTA file. Default is stdout.")
    
    parser.add_argument('-l', '--line_length',
                        action='store',
                        metavar='LENGTH',
                        type=int,
                        default=80,
                        help="Length of sequence lines in output FASTA file.")
    
    parser.add_argument('-v', '--verbose',
                        action='store_true',
                        help="Increase output verbosity.")
    
    return parser.parse_args()


def main() -> int:
    """Main function"""
    # Parse command line arguments
    args = parse_args()
    
    # Set verbosity
    if args.verbose:
        logger.setLevel(logging.DEBUG)
        logger.debug("Verbose output enabled")
    
    # Show input and output information if not using stdin/stdout
    if args.input_fastq is not sys.stdin:
        logger.info(f"Input file: {args.input_fastq.name}")
    if args.output_fasta is not sys.stdout:
        logger.info(f"Output file: {args.output_fasta.name}")
    
    # Perform the conversion
    seq_count = fastq_to_fasta(args.input_fastq, args.output_fasta, args.line_length)
    
    if seq_count >= 0:
        if args.verbose:
            logger.info(f"Successfully converted {seq_count} sequences")
        return 0
    else:
        logger.error("Conversion failed")
        return 1


if __name__ == '__main__':
    sys.exit(main())
    
