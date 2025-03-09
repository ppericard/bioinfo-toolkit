#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Test Data Generator

Utility functions for generating synthetic FASTQ and FASTA datasets for testing.
"""

import os
import random
import string
from pathlib import Path
from typing import Tuple, List, Dict, Optional

# Constants for generating random sequences
DNA_BASES = 'ACGT'

def generate_random_sequence(length: int) -> str:
    """Generate a random DNA sequence."""
    return ''.join(random.choice(DNA_BASES) for _ in range(length))

def generate_random_quality(length: int) -> str:
    """Generate a random quality string."""
    return ''.join(random.choice(string.ascii_letters) for _ in range(length))

def generate_read_id(read_number: int, pair_end: int = 0) -> str:
    """
    Generate a read ID in the format used by many sequencing platforms.
    
    Args:
        read_number: The read number
        pair_end: 1 for first in pair, 2 for second in pair, 0 for unpaired
        
    Returns:
        A formatted read ID string
    """
    if pair_end > 0:
        return f"@READ{read_number}/{pair_end}"
    else:
        return f"@READ{read_number}"

def write_fastq_record(file_handle, read_id: str, sequence: str, quality: str) -> None:
    """Write a FASTQ record to a file."""
    file_handle.write(f"{read_id}\n{sequence}\n+\n{quality}\n")

def create_paired_fastq_files(
    output_dir: Path,
    total_reads: int,
    paired_percent: float = 90.0,
    read_length: int = 150
) -> Tuple[Path, Path]:
    """
    Create synthetic paired FASTQ files with a specified percentage of paired reads.
    
    Args:
        output_dir: Directory to write output files
        total_reads: Total number of reads to generate
        paired_percent: Percentage of reads that should be paired (0-100)
        read_length: Length of reads to generate
        
    Returns:
        Tuple of (R1 file path, R2 file path)
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Calculate number of paired reads
    paired_reads = int(total_reads * (paired_percent / 100))
    unpaired_reads = total_reads - paired_reads
    
    # Create output files
    r1_file = output_dir / "reads_1.fastq"
    r2_file = output_dir / "reads_2.fastq"
    
    with open(r1_file, 'w') as r1, open(r2_file, 'w') as r2:
        # Generate paired reads
        for i in range(paired_reads):
            # Read 1
            read_id = generate_read_id(i, 1)
            sequence = generate_random_sequence(read_length)
            quality = generate_random_quality(read_length)
            write_fastq_record(r1, read_id, sequence, quality)
            
            # Read 2
            read_id = generate_read_id(i, 2)
            sequence = generate_random_sequence(read_length)
            quality = generate_random_quality(read_length)
            write_fastq_record(r2, read_id, sequence, quality)
        
        # Generate unpaired reads for R1
        unpaired_r1 = unpaired_reads // 2
        for i in range(paired_reads, paired_reads + unpaired_r1):
            read_id = generate_read_id(i, 1)
            sequence = generate_random_sequence(read_length)
            quality = generate_random_quality(read_length)
            write_fastq_record(r1, read_id, sequence, quality)
            
        # Generate unpaired reads for R2
        unpaired_r2 = unpaired_reads - unpaired_r1
        for i in range(paired_reads + unpaired_r1, total_reads):
            read_id = generate_read_id(i, 2)
            sequence = generate_random_sequence(read_length)
            quality = generate_random_quality(read_length)
            write_fastq_record(r2, read_id, sequence, quality)
    
    return r1_file, r2_file

def create_fasta_dataset(
    output_path: Path,
    num_sequences: int,
    min_length: int = 100,
    max_length: int = 1000
) -> Path:
    """
    Create a synthetic FASTA file with random sequences.
    
    Args:
        output_path: Path to write output file
        num_sequences: Number of sequences to generate
        min_length: Minimum sequence length
        max_length: Maximum sequence length
        
    Returns:
        Path to the created FASTA file
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_path.parent, exist_ok=True)
    
    with open(output_path, 'w') as f:
        for i in range(num_sequences):
            # Generate random sequence length
            seq_length = random.randint(min_length, max_length)
            
            # Generate sequence
            sequence = generate_random_sequence(seq_length)
            
            # Write to file
            f.write(f">SEQ{i} length={seq_length}\n")
            
            # Write sequence with line wrapping at 80 characters
            for j in range(0, len(sequence), 80):
                f.write(sequence[j:j+80] + "\n")
    
    return output_path 