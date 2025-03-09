#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Get Pairs

Description: Get separately paired reads and singletons 
             from two FASTQ files (left and right)

Examples:
  get_pairs.py -l file1.fastq -r file2.fastq -o output_dir
  get_pairs.py file1.fastq file2.fastq

-----------------------------------------------------------------------

Author: This software is written and maintained by Pierre Pericard
(pierre.pericard@ed.univ-lille1.fr)
Created: 2012-11-09
Modified: 2023-06-09
Licence: GNU GPL 3.0

Copyright 2013-2023 Pierre Pericard

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
from typing import Dict, Set, Optional

# Import local modules
try:
    # Try to import from the src package
    from bioinfotoolkit.utils.fastq_utils import build_read_id_set, process_fastq_file, extract_read_id, open_file
    from bioinfotoolkit.utils.bioinfo_logger import get_logger
except ImportError:
    # Fall back to relative imports for standalone usage
    sys.path.append(str(Path(__file__).parent.parent.parent.parent))
    from bioinfotoolkit.utils.fastq_utils import build_read_id_set, process_fastq_file, extract_read_id, open_file
    from bioinfotoolkit.utils.bioinfo_logger import get_logger


def get_pairs(
    left_file: str,
    right_file: str,
    output_dir: str,
    galaxy_mode: bool = False,
    verbose: bool = False
) -> Dict[str, Dict[str, int]]:
    """
    Process paired-end FASTQ files to separate paired and unpaired reads.
    
    Args:
        left_file: Path to the left reads FASTQ file
        right_file: Path to the right reads FASTQ file
        output_dir: Directory to write output files
        galaxy_mode: Whether to use Galaxy naming conventions
        verbose: Enable verbose output
        
    Returns:
        Dictionary with statistics for both files
    """
    # Setup logger
    logger = get_logger("get_pairs")
    log_level = logging.DEBUG if verbose else logging.INFO
    logger.set_level(log_level)
    
    # Create output directory if it doesn't exist
    output_dir_path = Path(output_dir)
    output_dir_path.mkdir(parents=True, exist_ok=True)
    
    # Log parameters
    logger.info(f"Left reads file: {left_file}")
    logger.info(f"Right reads file: {right_file}")
    logger.info(f"Output directory: {output_dir}")
    
    # Determine output filenames
    output_files = {}
    for file_path, suffix in [(left_file, "left"), (right_file, "right")]:
        if galaxy_mode:
            base_filename = suffix
            extension = "fastq"
        else:
            file_path_obj = Path(file_path)
            base_filename = file_path_obj.stem
            extension = file_path_obj.suffix.lstrip('.')
            if extension == 'gz':  # Handle .fastq.gz extensions
                extension = file_path_obj.suffixes[-2].lstrip('.') + '.gz'
                base_filename = file_path_obj.name[:-(len(extension)+1)]
        
        # Use string paths instead of Path objects for compatibility
        output_files[suffix] = {
            "paired": str(output_dir_path / f"{suffix}.paired.{extension}"),
            "unpaired": str(output_dir_path / f"{suffix}.unpaired.{extension}")
        }
    
    # STEP 1: Build an index of read IDs from the left file
    logger.start_progress(f"Building index from {left_file}")
    left_read_ids = build_read_id_set(
        left_file, 
        progress_callback=logger.get_progress_callback("Indexing left file")
    )
    logger.finish_progress()
    logger.info(f"Found {len(left_read_ids):,} unique read IDs in left file")
    
    # STEP 2: Process the right file using the left read IDs
    logger.start_progress(f"Processing {right_file}")
    right_stats = process_fastq_file(
        right_file,
        output_files["right"]["paired"],
        output_files["right"]["unpaired"],
        left_read_ids,
        invert_match=False,
        progress_callback=logger.get_progress_callback("Processing right file")
    )
    logger.finish_progress()
    
    # STEP 3: Extract paired read IDs from the processed right file
    # This ensures we only keep reads paired in both directions
    logger.start_progress("Identifying paired reads")
    
    # Optimization: Instead of reading the entire output file again,
    # we can use the intersection of left_read_ids and right_read_ids
    # that were actually found in the right file
    paired_read_ids = set()
    with open_file(output_files["right"]["paired"]) as f:
        record_count = 0
        for line_idx, line in enumerate(f):
            if line_idx % 4 == 0:  # Header line
                read_id = extract_read_id(line)
                paired_read_ids.add(read_id)
                record_count += 1
                
                # Periodically update progress
                if record_count % 10000 == 0:
                    logger.update_progress(record_count)
    
    logger.finish_progress()
    logger.info(f"Found {len(paired_read_ids):,} paired reads")
    
    # STEP 4: Process the left file to create paired and unpaired outputs
    logger.start_progress(f"Processing {left_file}")
    left_stats = process_fastq_file(
        left_file,
        output_files["left"]["paired"],
        output_files["left"]["unpaired"],
        paired_read_ids,
        invert_match=False,
        progress_callback=logger.get_progress_callback("Processing left file")
    )
    logger.finish_progress()
    
    # Log summary statistics
    logger.info(f"Left file: {left_stats['paired']} paired, {left_stats['unpaired']} unpaired")
    logger.info(f"Right file: {right_stats['paired']} paired, {right_stats['unpaired']} unpaired")
    
    return {
        "left": left_stats,
        "right": right_stats
    }


def parse_args() -> argparse.Namespace:
    """Parse and validate command line arguments."""
    parser = argparse.ArgumentParser(
        description='Get separately paired reads and singletons from two FASTQ files (left and right)',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Input files
    input_group = parser.add_argument_group('Input Files')
    input_group.add_argument('-l', '--left', 
                        dest='left_file',
                        help='Left reads FASTQ file')
    
    input_group.add_argument('-r', '--right', 
                        dest='right_file',
                        help='Right reads FASTQ file')
    
    # Optional arguments
    output_group = parser.add_argument_group('Output Options')
    output_group.add_argument('-o', '--output-dir', 
                        dest='output_dir',
                        default='.',
                        help='Output directory for paired and unpaired files')
    
    output_group.add_argument('--galaxy', 
                        action="store_true", 
                        default=False, 
                        help="Enable Galaxy mode (uses fixed output filenames)")
    
    # Performance options
    perf_group = parser.add_argument_group('Performance Options')
    perf_group.add_argument('--chunk-size',
                        type=int,
                        default=10000,
                        help='Number of records to process in each batch')
    
    # Logging options
    log_group = parser.add_argument_group('Logging Options')
    log_group.add_argument('-v', '--verbose',
                        action='store_true',
                        help='Enable verbose output')
    
    log_group.add_argument('-q', '--quiet',
                        action='store_true',
                        help='Suppress progress messages')
    
    # Legacy positional arguments
    parser.add_argument('leftreads', 
                        nargs='?', 
                        help=argparse.SUPPRESS)
    
    parser.add_argument('rightreads', 
                        nargs='?', 
                        help=argparse.SUPPRESS)
    
    args = parser.parse_args()
    
    # Handle legacy positional arguments
    if args.leftreads and args.rightreads:
        if not args.left_file:
            args.left_file = args.leftreads
        if not args.right_file:
            args.right_file = args.rightreads
    
    # Validate required arguments
    if not args.left_file or not args.right_file:
        parser.error("Must provide both left and right FASTQ files")
    
    return args


def main() -> int:
    """Main function."""
    # Parse command line arguments
    args = parse_args()
    
    try:
        # Process the files
        get_pairs(
            left_file=args.left_file,
            right_file=args.right_file,
            output_dir=args.output_dir,
            galaxy_mode=args.galaxy,
            verbose=args.verbose
        )
        return 0
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(main())
