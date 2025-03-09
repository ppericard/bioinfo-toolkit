#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
fasta_length_filter

Description: Filter a fasta file based on sequence length

Examples:
  fasta_length_filter.py -i input.fa -o output.fa -m 300
  fasta_length_filter.py -i input.fa -o output.fa -M 1000
  fasta_length_filter.py -i input.fa -o output.fa -m 300 -M 1000

-----------------------------------------------------------------------

Author: This software is written and maintained by Pierre Pericard
(pierre.pericard@ed.univ-lille1.fr)
Created: 2016-04-12
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
import logging
from typing import Generator, Tuple, TextIO, Optional, Dict
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


def read_fasta_file_handle(fasta_file_handle: TextIO) -> Generator[Tuple[str, str], None, None]:
    """
    Parse a fasta file and yield sequences as (header, sequence) tuples.
    
    Args:
        fasta_file_handle: An open file handle for the FASTA file
        
    Yields:
        Tuple containing (header, sequence)
    """
    # Variables initialization
    header = ''
    seqlines = []
    sequence_nb = 0
    
    try:
        # Reading input file
        for line in fasta_file_handle:
            line = line.strip()
            if not line:
                continue
                
            if line[0] == '>':
                # Yield the last read header and sequence
                if sequence_nb:
                    yield (header, ''.join(seqlines))
                    seqlines = []
                # Get header
                header = line[1:].rstrip()
                sequence_nb += 1
            else:
                # Concatenate sequence
                seqlines.append(line)
                
        # Yield the input file last sequence
        if header:
            yield (header, ''.join(seqlines))
            
    except Exception as e:
        logger.error(f"Error reading FASTA file: {e}")
        raise
    finally:
        # Don't close the file if it's stdin
        if fasta_file_handle is not sys.stdin:
            fasta_file_handle.close()


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


def filter_fasta_by_length(
    input_file: TextIO, 
    output_file: TextIO, 
    min_length: int = 0, 
    max_length: Optional[int] = None
) -> Dict[str, int]:
    """
    Filter sequences in a FASTA file based on their length.
    
    Args:
        input_file: Open file handle for input FASTA
        output_file: Open file handle for output FASTA
        min_length: Minimum sequence length to keep (default: 0)
        max_length: Maximum sequence length to keep (default: None, no maximum)
        
    Returns:
        Dictionary with statistics about the filtering process
    """
    stats = {
        'total': 0,
        'passed': 0,
        'filtered': 0
    }
    
    try:
        for header, sequence in read_fasta_file_handle(input_file):
            stats['total'] += 1
            seq_len = len(sequence)
            
            # Check if sequence passes length filters
            if seq_len < min_length:
                stats['filtered'] += 1
                continue
                
            if max_length is not None and seq_len > max_length:
                stats['filtered'] += 1
                continue
                
            # Write sequence to output file
            output_file.write(f">{header}\n{format_seq(sequence)}\n")
            stats['passed'] += 1
            
        return stats
        
    except Exception as e:
        logger.error(f"Error during filtering: {e}")
        return stats


def parse_args() -> argparse.Namespace:
    """Parse and return command line arguments"""
    parser = argparse.ArgumentParser(
        description='Filter a FASTA file based on sequence length.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('-i', '--input_fasta', 
                        metavar='INPUT', 
                        type=argparse.FileType('r'), 
                        default=sys.stdin,
                        help='Input FASTA file. Default is stdin.')
                        
    parser.add_argument('-o', '--output_fasta', 
                        metavar='OUTPUT', 
                        type=argparse.FileType('w'), 
                        default=sys.stdout,
                        help='Output FASTA file. Default is stdout.')
                        
    parser.add_argument('-m', '--min_length', 
                        metavar='MIN',
                        type=int, 
                        default=0,
                        help='Minimum sequence length to keep.')
                        
    parser.add_argument('-M', '--max_length', 
                        metavar='MAX',
                        type=int, 
                        default=None,
                        help='Maximum sequence length to keep. Default is no maximum.')
                        
    parser.add_argument('-v', '--verbose',
                        action='store_true',
                        help='Increase output verbosity.')
    
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
    if args.input_fasta is not sys.stdin:
        logger.info(f"Input file: {args.input_fasta.name}")
    if args.output_fasta is not sys.stdout:
        logger.info(f"Output file: {args.output_fasta.name}")
    
    # Log filter parameters
    logger.info(f"Minimum length: {args.min_length}")
    if args.max_length is not None:
        logger.info(f"Maximum length: {args.max_length}")
    else:
        logger.info("No maximum length specified")
    
    # Perform the filtering
    stats = filter_fasta_by_length(
        args.input_fasta, 
        args.output_fasta, 
        args.min_length, 
        args.max_length
    )
    
    # Log results
    if args.verbose:
        logger.info(f"Total sequences processed: {stats['total']}")
        logger.info(f"Sequences passed filter: {stats['passed']}")
        logger.info(f"Sequences filtered out: {stats['filtered']}")
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
