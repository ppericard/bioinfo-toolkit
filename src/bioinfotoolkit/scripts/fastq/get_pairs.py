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
from typing import Dict, Optional

from bioinfotoolkit.scripts.fastq.get_pairs_implementations import IMPLEMENTATIONS

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

def get_pairs(
    left_file: str,
    right_file: str,
    output_dir: str,
    implementation: str = 'v2',  # Default to v2 (improved in-memory implementation)
    compress: bool = False,
    verbose: bool = False,
    **kwargs
) -> Dict[str, Dict[str, int]]:
    """
    Get separately paired reads and singletons from two FASTQ files.
    
    Args:
        left_file: Path to the left (R1) FASTQ file
        right_file: Path to the right (R2) FASTQ file
        output_dir: Directory to write output files
        implementation: Which implementation to use ('v1', 'v2', or 'v3')
        compress: Whether to compress output files
        verbose: Whether to print verbose output
        **kwargs: Additional arguments to pass to the implementation
        
    Returns:
        Dictionary with counts of paired and singleton reads
    """
    # Convert paths to Path objects
    left_path = Path(left_file)
    right_path = Path(right_file)
    output_path = Path(output_dir)
    
    # Check if files exist
    if not left_path.exists():
        raise FileNotFoundError(f"Left file not found: {left_file}")
    
    if not right_path.exists():
        raise FileNotFoundError(f"Right file not found: {right_file}")
    
    # Check if implementation exists
    if implementation not in IMPLEMENTATIONS:
        raise ValueError(f"Unknown implementation: {implementation}. Available: {', '.join(IMPLEMENTATIONS.keys())}")
    
    # Get the implementation class
    impl_class = IMPLEMENTATIONS[implementation]
    
    # Run the implementation
    if verbose:
        logger.info(f"Using implementation: {implementation}")
    
    counts = impl_class.process(
        left_path,
        right_path,
        output_path,
        compress=compress,
        verbose=verbose,
        **kwargs
    )
    
    # Return the results
    return {
        'counts': counts,
        'files': {
            'paired_1': str(output_path / f"paired_1.fastq{'.gz' if compress else ''}"),
            'paired_2': str(output_path / f"paired_2.fastq{'.gz' if compress else ''}"),
            'singleton_1': str(output_path / f"singleton_1.fastq{'.gz' if compress else ''}"),
            'singleton_2': str(output_path / f"singleton_2.fastq{'.gz' if compress else ''}")
        }
    }

def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Get separately paired reads and singletons from two FASTQ files",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Input files
    input_group = parser.add_argument_group('Input')
    input_group.add_argument(
        '-l', '--left',
        dest='left_file',
        required=True,
        help='Left FASTQ file'
    )
    
    input_group.add_argument(
        '-r', '--right',
        dest='right_file',
        required=True,
        help='Right FASTQ file'
    )
    
    # Output options
    output_group = parser.add_argument_group('Output')
    output_group.add_argument(
        '-o', '--output-dir',
        dest='output_dir',
        default='get_pairs_output',
        help='Output directory'
    )
    
    output_group.add_argument(
        '-z', '--compress',
        action='store_true',
        help='Compress output files with gzip'
    )
    
    # Implementation options
    impl_group = parser.add_argument_group('Implementation')
    impl_group.add_argument(
        '-i', '--implementation',
        choices=list(IMPLEMENTATIONS.keys()),
        default='v2',
        help='Which implementation to use'
    )
    
    # V3-specific options
    v3_group = parser.add_argument_group('V3 Implementation Options')
    v3_group.add_argument(
        '--chunk-size',
        type=int,
        default=1000000,
        help='Number of reads to process at once (V3 only)'
    )
    
    v3_group.add_argument(
        '--temp-dir',
        help='Directory for temporary files (V3 only)'
    )
    
    # Other options
    other_group = parser.add_argument_group('Other')
    other_group.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Print verbose output'
    )
    
    # Galaxy compatibility
    other_group.add_argument(
        '-g', '--galaxy-mode',
        action='store_true',
        help='Galaxy mode (for compatibility)'
    )
    
    # Parse positional arguments if provided
    if len(sys.argv) == 3 and not sys.argv[1].startswith('-') and not sys.argv[2].startswith('-'):
        sys.argv = [sys.argv[0], '-l', sys.argv[1], '-r', sys.argv[2]]
    
    return parser.parse_args()

def main() -> int:
    """Main function."""
    args = parse_args()
    
    try:
        # Create output directory if it doesn't exist
        output_dir = Path(args.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Prepare kwargs for the implementation
        kwargs = {}
        if args.implementation == 'v3':
            kwargs['chunk_size'] = args.chunk_size
            if args.temp_dir:
                kwargs['temp_dir'] = Path(args.temp_dir)
        
        # Run get_pairs
        result = get_pairs(
            args.left_file,
            args.right_file,
            args.output_dir,
            implementation=args.implementation,
            compress=args.compress,
            verbose=args.verbose,
            **kwargs
        )
        
        # Print summary
        counts = result['counts']
        print("\nSummary:")
        print(f"  Paired reads: {counts['paired']}")
        print(f"  Singleton reads (left): {counts['singleton_1']}")
        print(f"  Singleton reads (right): {counts['singleton_2']}")
        print(f"  Total reads (left): {counts['total_1']}")
        print(f"  Total reads (right): {counts['total_2']}")
        
        # Print output files
        files = result['files']
        print("\nOutput files:")
        print(f"  Paired reads (left): {files['paired_1']}")
        print(f"  Paired reads (right): {files['paired_2']}")
        print(f"  Singleton reads (left): {files['singleton_1']}")
        print(f"  Singleton reads (right): {files['singleton_2']}")
        
        return 0
    
    except Exception as e:
        logger.error(f"Error: {e}")
        return 1

if __name__ == "__main__":
    sys.exit(main())
