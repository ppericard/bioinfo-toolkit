#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Get Pairs

Description: Get separately paired reads and singletons 
             from two FASTQ files (left and right)

Examples:
  get_pairs.py -l file1.fastq -r file2.fastq -o output_dir
  get_pairs.py file1.fastq file2.fastq

Available implementations:
  - v1: Original implementation from 2012-2016 (simple, minimal dependencies)
  - v2: Improved implementation using in-memory dictionaries
  - v3: Memory-optimized implementation using disk-based approach

-----------------------------------------------------------------------

Author: This software is written and maintained by Pierre Pericard
(pierre.pericard@ed.univ-lille1.fr)
Created: 2012-11-09
Modified: 2023-06-09
Licence: GNU GPL 3.0

Copyright 2012-2023 Pierre Pericard

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
from pathlib import Path
from typing import Dict, Union

# Import implementations
from bioinfotoolkit.scripts.fastq.get_pairs_v1 import GetPairsV1
from bioinfotoolkit.scripts.fastq.get_pairs_v2 import GetPairsV2
from bioinfotoolkit.scripts.fastq.get_pairs_v3 import GetPairsV3
from bioinfotoolkit.scripts.fastq.get_pairs_common import ensure_directory

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

# Implementation mapping
IMPLEMENTATIONS = {
    'v1': GetPairsV1,
    'v2': GetPairsV2,
    'v3': GetPairsV3
}

def get_pairs(
    left_file: Union[str, Path],
    right_file: Union[str, Path],
    output_dir: Union[str, Path],
    implementation: str = 'v2',  # Default to v2 (improved in-memory implementation)
    compress: bool = False,
    verbose: bool = False,
    **kwargs
) -> Dict[str, Dict[str, int]]:
    """
    Process paired-end FASTQ files to separate paired reads and singletons.
    
    Args:
        left_file: Path to the left reads FASTQ file
        right_file: Path to the right reads FASTQ file
        output_dir: Directory to write the output files
        implementation: Which implementation to use ('v1', 'v2', or 'v3')
        compress: Whether to compress the output files
        verbose: Whether to print verbose progress information
        **kwargs: Additional arguments to pass to the implementation
        
    Returns:
        Dictionary with statistics about the processing
    """
    # Convert to Path objects
    left_path = Path(left_file)
    right_path = Path(right_file)
    output_path = Path(output_dir)
    
    # Create output directory if it doesn't exist
    ensure_directory(output_path)
    
    # Get the implementation
    if implementation not in IMPLEMENTATIONS:
        logger.error(f"Unknown implementation: {implementation}")
        logger.error(f"Available implementations: {', '.join(IMPLEMENTATIONS.keys())}")
        return {}
    
    implementation_class = IMPLEMENTATIONS[implementation]
    logger.info(f"Using implementation: {implementation}")
    
    # Run the implementation
    stats = implementation_class.process(
        left_path, right_path, output_path, 
        compress=compress, verbose=verbose, **kwargs
    )
    
    return {'stats': stats}

def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Get separately paired reads and singletons from two FASTQ files (left and right)"
    )
    
    # Input files
    input_group = parser.add_argument_group('Input')
    input_files = input_group.add_mutually_exclusive_group(required=True)
    input_files.add_argument(
        'input_files', nargs='*', metavar='file', 
        help='Input FASTQ files (left and right)'
    )
    input_files.add_argument(
        '-l', '--left', dest='left_file',
        help='Left reads FASTQ file'
    )
    input_files.add_argument(
        '-r', '--right', dest='right_file',
        help='Right reads FASTQ file'
    )
    
    # Output options
    output_group = parser.add_argument_group('Output')
    output_group.add_argument(
        '-o', '--outdir', dest='output_dir', default='get_pairs_output',
        help='Output directory (default: get_pairs_output)'
    )
    output_group.add_argument(
        '-z', '--gzip', dest='compress', action='store_true',
        help='Compress output files with gzip'
    )
    
    # Implementation selection
    impl_group = parser.add_argument_group('Implementation')
    impl_group.add_argument(
        '-i', '--implementation', dest='implementation', default='v2',
        choices=IMPLEMENTATIONS.keys(),
        help='Implementation to use (default: v2)'
    )
    
    # Logging options
    log_group = parser.add_argument_group('Logging')
    log_group.add_argument(
        '-v', '--verbose', dest='verbose', action='store_true',
        help='Print verbose progress information'
    )
    
    args = parser.parse_args()
    
    # Handle positional arguments if provided
    if args.input_files and len(args.input_files) >= 2:
        args.left_file = args.input_files[0]
        args.right_file = args.input_files[1]
    
    # Validate input files
    if not args.left_file or not args.right_file:
        parser.error('Both left and right FASTQ files are required')
    
    return args

def main() -> int:
    """Main entry point for the script."""
    args = parse_args()
    
    try:
        results = get_pairs(
            args.left_file,
            args.right_file,
            args.output_dir,
            implementation=args.implementation,
            compress=args.compress,
            verbose=args.verbose
        )
        
        if not results:
            return 1
        
        stats = results['stats']
        logger.info("Processing complete!")
        logger.info(f"Paired reads: {stats['paired']}")
        logger.info(f"Left singletons: {stats['left_singletons']}")
        logger.info(f"Right singletons: {stats['right_singletons']}")
        
        return 0
    
    except Exception as e:
        logger.error(f"Error: {e}")
        return 1

if __name__ == "__main__":
    sys.exit(main())
