#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Get Pairs (Memory-Optimized Version)

Description: Get separately paired reads and singletons from two FASTQ files (left and right)
            using a disk-based approach that minimizes memory usage.

Examples:
  get_pairs_v3.py -l file1.fastq -r file2.fastq -o output_dir
  get_pairs_v3.py file1.fastq file2.fastq

-----------------------------------------------------------------------

Author: This software is written and maintained by Pierre Pericard
(pierre.pericard@ed.univ-lille1.fr)
Created: 2023-06-09
Licence: GNU GPL 3.0

Copyright 2023 Pierre Pericard

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
import tempfile
import gzip
import shutil
import time
from pathlib import Path
from typing import Dict, List, Set, Tuple, Optional, TextIO, Union, Generator, BinaryIO

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger("get_pairs_v3")


def open_file(filename: Union[str, Path], mode: str = 'r') -> TextIO:
    """
    Open a file for reading or writing, handling gzipped files automatically.
    
    Args:
        filename: Path to the file
        mode: File mode ('r' for reading, 'w' for writing)
        
    Returns:
        File handle
    """
    # Convert to Path object
    file_path = Path(filename)
    
    # Check if the file is gzipped
    if str(file_path).endswith('.gz'):
        return gzip.open(file_path, mode + 't')
    else:
        return open(file_path, mode)


def extract_read_id(header_line: str) -> str:
    """
    Extract the read ID from a FASTQ header line.
    
    Args:
        header_line: The FASTQ header line starting with '@'
        
    Returns:
        Read ID without the '@' prefix and direction suffix (/1 or /2)
    """
    # Remove leading '@' and strip whitespace
    header = header_line.strip()
    if header.startswith('@'):
        header = header[1:]
    
    # Get first part of the header (before whitespace)
    header = header.split()[0]
    
    # Remove /1 or /2 suffix if present
    if '/' in header:
        header = header.split('/')[0]
    
    return header


def extract_and_sort_ids(fastq_file: str, temp_dir: str) -> str:
    """
    Extract read IDs from a FASTQ file, sort them, and write to a temporary file.
    This is memory-efficient as it only stores one read ID at a time.
    
    Args:
        fastq_file: Path to the FASTQ file
        temp_dir: Directory to write temporary files
        
    Returns:
        Path to the sorted IDs file
    """
    start_time = time.time()
    read_count = 0
    
    # Create a temporary file for unsorted IDs
    unsorted_ids_file = Path(temp_dir) / f"unsorted_ids_{Path(fastq_file).stem}.txt"
    
    logger.info(f"Extracting read IDs from {fastq_file}")
    
    # Extract read IDs to unsorted temporary file
    with open_file(fastq_file) as fastq, open(unsorted_ids_file, 'w') as ids_file:
        line_number = 0
        for line in fastq:
            line_number += 1
            if line_number % 4 == 1:  # Header lines
                read_id = extract_read_id(line)
                ids_file.write(f"{read_id}\n")
                read_count += 1
                
                # Print progress every 100,000 reads
                if read_count % 100000 == 0:
                    elapsed = time.time() - start_time
                    logger.info(f"Processed {read_count:,} reads ({read_count/elapsed:.2f} reads/sec)")
    
    # Sort the IDs file externally (memory efficient)
    logger.info(f"Sorting {read_count:,} read IDs...")
    sorted_ids_file = Path(temp_dir) / f"sorted_ids_{Path(fastq_file).stem}.txt"
    
    # Use Unix sort if available (faster and more memory efficient)
    # On Windows, fallback to Python sorting
    try:
        import subprocess
        # Check if sort command exists
        try:
            subprocess.run(["sort", "--version"], capture_output=True, check=True)
            # Use external sort command
            subprocess.run(["sort", "-o", str(sorted_ids_file), str(unsorted_ids_file)], check=True)
        except (subprocess.SubprocessError, FileNotFoundError):
            # Fallback to Python sorting
            raise ImportError("External sort command not available")
    except ImportError:
        # Python-based sorting using minimal memory
        logger.info("Using Python-based sorting (may be slower)")
        
        # Read in small chunks, sort, and write to temporary files
        chunk_size = 100000  # Adjust based on available memory
        temp_files = []
        
        with open(unsorted_ids_file, 'r') as in_file:
            chunk = []
            chunk_num = 0
            
            for line in in_file:
                read_id = line.strip()
                chunk.append(read_id)
                
                if len(chunk) >= chunk_size:
                    chunk_num += 1
                    chunk_file = Path(temp_dir) / f"chunk_{chunk_num}.txt"
                    temp_files.append(chunk_file)
                    
                    # Sort this chunk and write to temporary file
                    chunk.sort()
                    with open(chunk_file, 'w') as out_file:
                        for read_id in chunk:
                            out_file.write(f"{read_id}\n")
                    
                    chunk = []
            
            # Handle the last chunk
            if chunk:
                chunk_num += 1
                chunk_file = Path(temp_dir) / f"chunk_{chunk_num}.txt"
                temp_files.append(chunk_file)
                
                chunk.sort()
                with open(chunk_file, 'w') as out_file:
                    for read_id in chunk:
                        out_file.write(f"{read_id}\n")
        
        # Merge sorted temporary files
        if temp_files:
            with open(sorted_ids_file, 'w') as out_file:
                file_handles = [open(f, 'r') for f in temp_files]
                lines = [f.readline().strip() for f in file_handles]
                
                # Merge until all files are exhausted
                while any(lines):
                    # Find the minimum value among the current lines
                    min_idx = -1
                    min_val = None
                    
                    for i, val in enumerate(lines):
                        if val and (min_val is None or val < min_val):
                            min_idx = i
                            min_val = val
                    
                    if min_idx >= 0:
                        # Write the minimum value
                        out_file.write(f"{min_val}\n")
                        
                        # Get the next line from the corresponding file
                        lines[min_idx] = file_handles[min_idx].readline().strip()
                
                # Close all file handles
                for fh in file_handles:
                    fh.close()
                
                # Delete temporary chunk files
                for f in temp_files:
                    os.unlink(f)
    
    # Remove the unsorted file
    os.unlink(unsorted_ids_file)
    
    elapsed = time.time() - start_time
    logger.info(f"Extracted and sorted {read_count:,} read IDs in {elapsed:.2f} seconds")
    
    return str(sorted_ids_file)


def merge_id_files(left_ids_file: str, right_ids_file: str, temp_dir: str) -> Tuple[int, int, int]:
    """
    Merge two sorted ID files to identify paired and unpaired reads.
    This is memory-efficient as it only reads one line from each file at a time.
    
    Args:
        left_ids_file: Path to the sorted IDs from the left file
        right_ids_file: Path to the sorted IDs from the right file
        temp_dir: Directory to write output files
        
    Returns:
        Tuple of (path to paired IDs file, path to left-only IDs file, path to right-only IDs file)
    """
    start_time = time.time()
    
    # Create output files
    paired_ids_file = Path(temp_dir) / "paired_ids.txt"
    left_only_ids_file = Path(temp_dir) / "left_only_ids.txt"
    right_only_ids_file = Path(temp_dir) / "right_only_ids.txt"
    
    # Counts for statistics
    paired_count = 0
    left_only_count = 0
    right_only_count = 0
    
    logger.info("Merging ID files to identify paired and unpaired reads...")
    
    with open(left_ids_file, 'r') as left_file, \
         open(right_ids_file, 'r') as right_file, \
         open(paired_ids_file, 'w') as paired_file, \
         open(left_only_ids_file, 'w') as left_only_file, \
         open(right_only_ids_file, 'w') as right_only_file:
        
        # Get initial lines
        left_line = left_file.readline().strip()
        right_line = right_file.readline().strip()
        
        # Process until both files are exhausted
        while left_line or right_line:
            # Both files have more lines
            if left_line and right_line:
                if left_line == right_line:
                    # IDs match - paired read
                    paired_file.write(f"{left_line}\n")
                    paired_count += 1
                    
                    # Move to next lines in both files
                    left_line = left_file.readline().strip()
                    right_line = right_file.readline().strip()
                elif left_line < right_line:
                    # Left ID is smaller - left-only read
                    left_only_file.write(f"{left_line}\n")
                    left_only_count += 1
                    
                    # Move to next line in left file
                    left_line = left_file.readline().strip()
                else:
                    # Right ID is smaller - right-only read
                    right_only_file.write(f"{right_line}\n")
                    right_only_count += 1
                    
                    # Move to next line in right file
                    right_line = right_file.readline().strip()
            
            # Only left file has more lines
            elif left_line:
                left_only_file.write(f"{left_line}\n")
                left_only_count += 1
                left_line = left_file.readline().strip()
            
            # Only right file has more lines
            elif right_line:
                right_only_file.write(f"{right_line}\n")
                right_only_count += 1
                right_line = right_file.readline().strip()
            
            # Progress logging (every million reads)
            total_processed = paired_count + left_only_count + right_only_count
            if total_processed % 1000000 == 0:
                elapsed = time.time() - start_time
                logger.info(f"Processed {total_processed:,} reads ({total_processed/elapsed:.2f} reads/sec)")
    
    elapsed = time.time() - start_time
    total_reads = paired_count + left_only_count + right_only_count
    
    logger.info(f"Merged ID files in {elapsed:.2f} seconds:")
    logger.info(f"  Paired reads: {paired_count:,}")
    logger.info(f"  Left-only reads: {left_only_count:,}")
    logger.info(f"  Right-only reads: {right_only_count:,}")
    
    return str(paired_ids_file), str(left_only_ids_file), str(right_only_ids_file)


def get_output_paths(file_path: str, output_dir: str, galaxy_mode: bool, suffix: str) -> Tuple[str, str]:
    """
    Determine output file paths based on input file, output directory, and mode.
    
    Args:
        file_path: Input file path
        output_dir: Output directory
        galaxy_mode: Whether to use Galaxy naming conventions
        suffix: Suffix for this file (left or right)
        
    Returns:
        Tuple of (paired output path, unpaired output path)
    """
    output_dir_path = Path(output_dir)
    
    if galaxy_mode:
        base_filename = suffix
        extension = "fastq"
    else:
        file_path_obj = Path(file_path)
        base_filename = file_path_obj.stem
        extension = file_path_obj.suffix.lstrip('.')
        
        # Handle .fastq.gz extensions
        if extension == 'gz':
            extension = file_path_obj.suffixes[-2].lstrip('.') + '.gz'
            base_filename = file_path_obj.name[:-(len(extension)+1)]
    
    paired_path = output_dir_path / f"{base_filename}.paired.{extension}"
    unpaired_path = output_dir_path / f"{base_filename}.unpaired.{extension}"
    
    return str(paired_path), str(unpaired_path)


def categorize_and_write_reads(
    fastq_file: str,
    paired_ids_file: str,
    unpaired_ids_file: str,
    paired_output: str,
    unpaired_output: str
) -> Dict[str, int]:
    """
    Categorize reads from a FASTQ file and write to paired and unpaired output files.
    Memory-efficient as it only loads the IDs of one type (paired or unpaired) at a time.
    
    Args:
        fastq_file: Path to the FASTQ file
        paired_ids_file: Path to file containing paired read IDs
        unpaired_ids_file: Path to file containing unpaired read IDs
        paired_output: Path to write paired reads
        unpaired_output: Path to write unpaired reads
        
    Returns:
        Dictionary with statistics about processed reads
    """
    start_time = time.time()
    
    # Load either paired or unpaired IDs depending on which set is smaller
    paired_size = os.path.getsize(paired_ids_file)
    unpaired_size = os.path.getsize(unpaired_ids_file)
    
    # Choose the smaller file to load into memory
    if paired_size <= unpaired_size:
        logger.info(f"Loading paired IDs into memory (smaller set)")
        with open(paired_ids_file, 'r') as f:
            ids = set(line.strip() for line in f)
        is_paired = lambda read_id: read_id in ids
    else:
        logger.info(f"Loading unpaired IDs into memory (smaller set)")
        with open(unpaired_ids_file, 'r') as f:
            ids = set(line.strip() for line in f)
        is_paired = lambda read_id: read_id not in ids
    
    logger.info(f"Loaded {len(ids):,} read IDs into memory")
    
    # Process the FASTQ file and write to output files
    stats = {
        'total': 0,
        'paired': 0,
        'unpaired': 0
    }
    
    logger.info(f"Processing {fastq_file} and writing to output files...")
    
    with open_file(fastq_file) as input_file, \
         open_file(paired_output, 'w') as paired_file, \
         open_file(unpaired_output, 'w') as unpaired_file:
        
        line_number = 0
        current_id = None
        current_record = []
        record_is_paired = False
        
        for line in input_file:
            line = line.strip()
            if not line:
                continue
            
            line_number += 1
            line_position = (line_number - 1) % 4
            
            if line_position == 0:
                # New record, process the previous one if it exists
                if current_record:
                    output_file = paired_file if record_is_paired else unpaired_file
                    output_file.write('\n'.join(current_record) + '\n')
                    
                    stats['total'] += 1
                    if record_is_paired:
                        stats['paired'] += 1
                    else:
                        stats['unpaired'] += 1
                    
                    if stats['total'] % 100000 == 0:
                        elapsed = time.time() - start_time
                        rate = stats['total'] / elapsed if elapsed > 0 else 0
                        logger.info(f"Processed {stats['total']:,} reads ({rate:.2f} reads/sec)")
                
                # Start new record
                current_id = extract_read_id(line)
                record_is_paired = is_paired(current_id)
                current_record = [line]
            else:
                # Continue current record
                current_record.append(line)
        
        # Process the last record
        if current_record:
            output_file = paired_file if record_is_paired else unpaired_file
            output_file.write('\n'.join(current_record) + '\n')
            
            stats['total'] += 1
            if record_is_paired:
                stats['paired'] += 1
            else:
                stats['unpaired'] += 1
    
    elapsed = time.time() - start_time
    rate = stats['total'] / elapsed if elapsed > 0 else 0
    
    logger.info(f"Processed {stats['total']:,} reads in {elapsed:.2f} seconds ({rate:.2f} reads/sec)")
    logger.info(f"  Paired: {stats['paired']:,}")
    logger.info(f"  Unpaired: {stats['unpaired']:,}")
    
    return stats


def get_pairs(
    left_file: str,
    right_file: str,
    output_dir: str,
    galaxy_mode: bool = False,
    keep_temp: bool = False,
    temp_dir: Optional[str] = None,
    verbose: bool = False
) -> Dict[str, Dict[str, int]]:
    """
    Process paired-end FASTQ files to separate paired and unpaired reads.
    Uses a disk-based approach to minimize memory usage.
    
    Args:
        left_file: Path to the left reads FASTQ file
        right_file: Path to the right reads FASTQ file
        output_dir: Directory to write output files
        galaxy_mode: Whether to use Galaxy naming conventions
        keep_temp: Whether to keep temporary files
        temp_dir: Directory for temporary files, or None to use system default
        verbose: Enable verbose output
        
    Returns:
        Dictionary with statistics for both files
    """
    # Set logging level
    if verbose:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
    
    # Create output directory
    output_dir_path = Path(output_dir)
    output_dir_path.mkdir(parents=True, exist_ok=True)
    
    # Create temporary directory
    using_temp_dir = temp_dir is None
    if using_temp_dir:
        temp_dir = tempfile.mkdtemp(prefix="get_pairs_")
    else:
        temp_dir = Path(temp_dir)
        os.makedirs(temp_dir, exist_ok=True)
    
    logger.info(f"Left reads file: {left_file}")
    logger.info(f"Right reads file: {right_file}")
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"Temporary directory: {temp_dir}")
    
    try:
        # Phase 1: Extract and sort read IDs from both files
        left_ids_file = extract_and_sort_ids(left_file, temp_dir)
        right_ids_file = extract_and_sort_ids(right_file, temp_dir)
        
        # Phase 2: Merge sorted ID files to identify paired and unpaired reads
        paired_ids_file, left_only_ids_file, right_only_ids_file = merge_id_files(
            left_ids_file, right_ids_file, temp_dir)
        
        # Get output file paths
        left_paired_output, left_unpaired_output = get_output_paths(
            left_file, output_dir, galaxy_mode, "left")
        
        right_paired_output, right_unpaired_output = get_output_paths(
            right_file, output_dir, galaxy_mode, "right")
        
        # Phase 3: Categorize and write reads
        logger.info("Processing left file...")
        left_stats = categorize_and_write_reads(
            left_file, paired_ids_file, left_only_ids_file, 
            left_paired_output, left_unpaired_output)
        
        logger.info("Processing right file...")
        right_stats = categorize_and_write_reads(
            right_file, paired_ids_file, right_only_ids_file,
            right_paired_output, right_unpaired_output)
        
        # Log results
        logger.info(f"Summary:")
        logger.info(f"  Left file: {left_stats['paired']:,} paired, {left_stats['unpaired']:,} unpaired")
        logger.info(f"  Right file: {right_stats['paired']:,} paired, {right_stats['unpaired']:,} unpaired")
        
        return {
            "left": left_stats,
            "right": right_stats
        }
    
    finally:
        # Clean up temporary files
        if not keep_temp and using_temp_dir:
            logger.debug(f"Cleaning up temporary directory: {temp_dir}")
            shutil.rmtree(temp_dir, ignore_errors=True)


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
    
    # Output options
    output_group = parser.add_argument_group('Output Options')
    output_group.add_argument('-o', '--output-dir', 
                        dest='output_dir',
                        default='.',
                        help='Output directory for paired and unpaired files')
    
    output_group.add_argument('--galaxy', 
                        action="store_true", 
                        default=False, 
                        help="Enable Galaxy mode (uses fixed output filenames)")
    
    # Temporary file options
    temp_group = parser.add_argument_group('Temporary Files')
    temp_group.add_argument('--temp-dir',
                        dest='temp_dir',
                        help='Directory for temporary files (default: system temp directory)')
    
    temp_group.add_argument('--keep-temp',
                        action='store_true',
                        help='Keep temporary files after completion')
    
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
    
    # Adjust logging level
    if args.quiet:
        logger.setLevel(logging.WARNING)
    
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
            keep_temp=args.keep_temp,
            temp_dir=args.temp_dir,
            verbose=args.verbose
        )
        return 0
    except Exception as e:
        logger.error(f"Error: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(main()) 