#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
FASTQ Utils

Description: Common utilities for processing FASTQ files
"""

import gzip
import os
from pathlib import Path
from typing import Generator, List, TextIO, Tuple, Union, Optional, Dict, Set


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
    Extract the read ID from a FASTQ header line, removing direction indicators.
    
    Args:
        header_line: The FASTQ header line starting with '@'
        
    Returns:
        Cleaned read ID (without '@' prefix or direction suffix)
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


def read_fastq_records(file_path: Union[str, Path]) -> Generator[Tuple[str, List[str]], None, None]:
    """
    Read a FASTQ file and yield records with their IDs.
    
    Args:
        file_path: Path to the FASTQ file
        
    Yields:
        Tuple of (read_id, [header, sequence, plus_line, quality])
    """
    with open_file(file_path) as f:
        record = []
        read_id = None
        record_count = 0
        
        for line_num, line in enumerate(f):
            line = line.strip()
            if not line:
                continue
                
            line_position = line_num % 4
            record.append(line)
            
            # The first line of each record is the header
            if line_position == 0:
                read_id = extract_read_id(line)
            
            # When we have all 4 lines of a record
            if line_position == 3:
                record_count += 1
                yield read_id, record
                record = []


def chunk_iterator(iterator, chunk_size: int = 1000):
    """
    Process an iterator in chunks for more efficient batch processing.
    
    Args:
        iterator: The source iterator to chunk
        chunk_size: Number of items to yield in each chunk
        
    Yields:
        Lists of items from the iterator, up to chunk_size in length
    """
    chunk = []
    for item in iterator:
        chunk.append(item)
        if len(chunk) >= chunk_size:
            yield chunk
            chunk = []
            
    # Yield the final chunk if not empty
    if chunk:
        yield chunk


def build_read_id_set(file_path: Union[str, Path], progress_callback=None) -> Set[str]:
    """
    Build a set of read IDs from a FASTQ file efficiently.
    
    Args:
        file_path: Path to the FASTQ file
        progress_callback: Optional callback function to report progress
        
    Returns:
        Set of read IDs
    """
    read_ids = set()
    record_count = 0
    
    # Process records in chunks for better performance
    for chunk_idx, chunk in enumerate(chunk_iterator(read_fastq_records(file_path), 10000)):
        for read_id, _ in chunk:
            read_ids.add(read_id)
            record_count += 1
        
        # Report progress if callback provided
        if progress_callback and chunk_idx % 10 == 0:
            progress_callback(record_count)
    
    return read_ids


def process_fastq_file(
    file_path: Union[str, Path],
    paired_output: Union[str, Path],
    unpaired_output: Union[str, Path],
    read_id_set: Set[str],
    invert_match: bool = False,
    progress_callback=None
) -> Dict[str, int]:
    """
    Process a FASTQ file and separate reads into paired and unpaired files.
    
    Args:
        file_path: Path to the input FASTQ file
        paired_output: Path to the output file for paired reads
        unpaired_output: Path to the output file for unpaired reads
        read_id_set: Set of read IDs to check against
        invert_match: If True, invert the match condition
        progress_callback: Optional callback function to report progress
        
    Returns:
        Dictionary with statistics about processed reads
    """
    stats = {
        'total': 0,
        'paired': 0,
        'unpaired': 0
    }
    
    with open_file(paired_output, 'w') as paired_file, open_file(unpaired_output, 'w') as unpaired_file:
        record_count = 0
        
        # Read and process records in chunks for better performance
        for chunk_idx, chunk in enumerate(chunk_iterator(read_fastq_records(file_path), 10000)):
            for read_id, record in chunk:
                # Determine if this read is paired
                is_in_set = read_id in read_id_set
                is_paired = is_in_set if not invert_match else not is_in_set
                
                # Write to appropriate file
                output_file = paired_file if is_paired else unpaired_file
                output_file.write('\n'.join(record) + '\n')
                
                # Update statistics
                stats['total'] += 1
                if is_paired:
                    stats['paired'] += 1
                else:
                    stats['unpaired'] += 1
                
                record_count += 1
            
            # Report progress if callback provided
            if progress_callback and chunk_idx % 10 == 0:
                progress_callback(record_count, stats)
    
    return stats 