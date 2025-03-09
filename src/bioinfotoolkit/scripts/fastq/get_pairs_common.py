#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Get Pairs Common Utilities

This module contains common utility functions used by different implementations
of the get_pairs algorithm for processing paired-end FASTQ files.
"""

import os
import logging
from pathlib import Path
from typing import Dict, TextIO, Set, Tuple, Optional, List, Generator, Union

from bioinfotoolkit.utils.fastq_utils import open_file, extract_read_id

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

# Common utility functions
def ensure_directory(directory: Path) -> None:
    """Ensure a directory exists, creating it if necessary."""
    directory.mkdir(parents=True, exist_ok=True)

def count_reads(fastq_file: Path) -> int:
    """
    Count the number of reads in a FASTQ file.
    
    Args:
        fastq_file: Path to the FASTQ file
        
    Returns:
        Number of reads in the file
    """
    count = 0
    with open_file(fastq_file) as f:
        for line in f:
            if line.startswith('@'):
                count += 1
                # Skip the next 3 lines (sequence, +, quality)
                next(f)
                next(f)
                next(f)
    return count

def write_fastq_record(file_handle: TextIO, header: str, sequence: str, plus_line: str, quality: str) -> None:
    """
    Write a FASTQ record to a file.
    
    Args:
        file_handle: File handle to write to
        header: FASTQ header line (starting with @)
        sequence: Sequence line
        plus_line: Plus line (usually just +)
        quality: Quality line
    """
    file_handle.write(f"{header}\n{sequence}\n{plus_line}\n{quality}\n")

def get_output_paths(output_dir: Path, compress: bool = False) -> Dict[str, Path]:
    """
    Get the paths to the output files.
    
    Args:
        output_dir: Output directory
        compress: Whether to compress the output files
        
    Returns:
        Dictionary with paths to the output files
    """
    extension = ".fastq.gz" if compress else ".fastq"
    
    return {
        'paired_1': output_dir / f"paired_1{extension}",
        'paired_2': output_dir / f"paired_2{extension}",
        'singleton_1': output_dir / f"singleton_1{extension}",
        'singleton_2': output_dir / f"singleton_2{extension}"
    }

def build_read_id_set_from_file(file_path: Path, verbose: bool = False) -> Set[str]:
    """
    Build a set of read IDs from a FASTQ file.
    
    Args:
        file_path: Path to the FASTQ file
        verbose: Whether to print verbose output
        
    Returns:
        Set of read IDs
    """
    if verbose:
        logger.info(f"Building read ID set from {file_path}")
    
    read_ids = set()
    
    with open_file(file_path) as f:
        line_count = 0
        for line in f:
            if line_count % 4 == 0:  # Header line
                read_id = extract_read_id(line)
                read_ids.add(read_id)
            line_count += 1
    
    if verbose:
        logger.info(f"Built read ID set with {len(read_ids)} unique IDs")
    
    return read_ids

def process_fastq_pairs(
    reads_1: Dict[str, List[str]],
    reads_2: Dict[str, List[str]],
    output_paths: Dict[str, Path],
    compress: bool = False,
    verbose: bool = False
) -> Dict[str, int]:
    """
    Process paired reads and write them to output files.
    
    Args:
        reads_1: Dictionary of read ID -> FASTQ record lines for left reads
        reads_2: Dictionary of read ID -> FASTQ record lines for right reads
        output_paths: Dictionary with paths to the output files
        compress: Whether to compress the output files
        verbose: Whether to print verbose output
        
    Returns:
        Dictionary with counts of paired and singleton reads
    """
    # Open output files
    with open_file(output_paths['paired_1'], 'w') as paired_1_file, \
         open_file(output_paths['paired_2'], 'w') as paired_2_file, \
         open_file(output_paths['singleton_1'], 'w') as singleton_1_file, \
         open_file(output_paths['singleton_2'], 'w') as singleton_2_file:
        
        # Process paired reads
        paired_count = 0
        singleton_1_count = 0
        singleton_2_count = 0
        
        # Find common reads
        common_ids = set(reads_1.keys()) & set(reads_2.keys())
        
        if verbose:
            logger.info(f"Found {len(common_ids)} paired reads")
        
        # Write paired reads
        for read_id in common_ids:
            record_1 = reads_1[read_id]
            record_2 = reads_2[read_id]
            
            paired_1_file.write('\n'.join(record_1) + '\n')
            paired_2_file.write('\n'.join(record_2) + '\n')
            paired_count += 1
        
        # Write singletons from reads_1
        for read_id in set(reads_1.keys()) - common_ids:
            record = reads_1[read_id]
            singleton_1_file.write('\n'.join(record) + '\n')
            singleton_1_count += 1
        
        # Write singletons from reads_2
        for read_id in set(reads_2.keys()) - common_ids:
            record = reads_2[read_id]
            singleton_2_file.write('\n'.join(record) + '\n')
            singleton_2_count += 1
    
    return {
        'paired': paired_count,
        'left_singletons': singleton_1_count,
        'right_singletons': singleton_2_count,
        'total_left': paired_count + singleton_1_count,
        'total_right': paired_count + singleton_2_count
    }

def read_fastq_records_grouped(
    file_path: Path,
    chunk_size: int = 10000,
    verbose: bool = False
) -> Generator[Dict[str, List[str]], None, None]:
    """
    Read FASTQ records in chunks and group them by read ID.
    
    Args:
        file_path: Path to the FASTQ file
        chunk_size: Number of reads to process at once
        verbose: Whether to print verbose output
        
    Yields:
        Dictionary of read ID -> FASTQ record lines for each chunk
    """
    if verbose:
        logger.info(f"Reading records from {file_path} (chunk size: {chunk_size})")
    
    reads_chunk = {}
    read_count = 0
    
    with open_file(file_path) as f:
        current_record = []
        current_read_id = None
        
        for line_num, line in enumerate(f):
            line = line.rstrip()
            
            if line_num % 4 == 0:  # Header line
                if current_read_id and current_record:
                    reads_chunk[current_read_id] = current_record
                    read_count += 1
                    
                    if read_count >= chunk_size:
                        if verbose:
                            logger.info(f"Processed {read_count} reads")
                        yield reads_chunk
                        reads_chunk = {}
                        read_count = 0
                
                current_read_id = extract_read_id(line)
                current_record = [line]
            else:
                current_record.append(line)
        
        # Add the last record
        if current_read_id and current_record:
            reads_chunk[current_read_id] = current_record
            read_count += 1
    
    # Yield the last chunk
    if reads_chunk:
        if verbose:
            logger.info(f"Processed {read_count} reads (final chunk)")
        yield reads_chunk 