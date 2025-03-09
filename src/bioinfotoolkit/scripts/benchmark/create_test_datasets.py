#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Create Test Datasets

Description: Generate synthetic FASTQ datasets for testing get_pairs.py
             with different sizes and proportions of paired reads.
"""

import os
import random
import string
import argparse
import sys
from pathlib import Path
from typing import Tuple, List, Dict, Optional

# Constants for generating random sequences
DNA_BASES = 'ACGT'
QUAL_CHARS = ''.join(chr(x) for x in range(33, 74))  # ASCII 33-73 for quality scores


def generate_random_sequence(length: int) -> str:
    """Generate a random DNA sequence of specified length."""
    return ''.join(random.choice(DNA_BASES) for _ in range(length))


def generate_random_quality(length: int) -> str:
    """Generate a random quality string of specified length."""
    return ''.join(random.choice(QUAL_CHARS) for _ in range(length))


def generate_read_id(read_number: int, pair_end: int) -> str:
    """
    Generate a read ID in a format typical for paired-end sequencing.
    
    Args:
        read_number: The read number (unique identifier)
        pair_end: 1 for left read, 2 for right read
        
    Returns:
        Formatted read ID string
    """
    return f"@READID_{read_number:06d}/FLOWCELL:LANE:{read_number:06d}:{pair_end}"


def write_fastq_record(file_handle, read_id: str, sequence: str, quality: str) -> None:
    """Write a single FASTQ record to a file."""
    file_handle.write(f"{read_id}\n")
    file_handle.write(f"{sequence}\n")
    file_handle.write("+\n")
    file_handle.write(f"{quality}\n")


def create_paired_fastq_files(
    output_dir: Path,
    total_reads: int,
    paired_percent: float = 90.0,
    read_length: int = 150
) -> Tuple[Path, Path]:
    """
    Create paired FASTQ files with a specified percentage of paired reads.
    
    Args:
        output_dir: Directory to write output files
        total_reads: Total number of reads to generate (per file)
        paired_percent: Percentage of reads that should be paired (0-100)
        read_length: Length of each read
        
    Returns:
        Tuple of paths to the left and right FASTQ files
    """
    # Create output directory if it doesn't exist
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Calculate number of paired reads
    num_paired = int((total_reads * paired_percent) / 100)
    
    # File paths
    left_file_path = output_dir / f"reads_1_{total_reads}_{int(paired_percent)}.fastq"
    right_file_path = output_dir / f"reads_2_{total_reads}_{int(paired_percent)}.fastq"
    
    # Create the paired read IDs first
    paired_ids = list(range(1, num_paired + 1))
    
    # Generate different read IDs for the unpaired reads
    left_unpaired_ids = list(range(num_paired + 1, total_reads + 1))
    right_unpaired_ids = list(range(total_reads + 1, total_reads + (total_reads - num_paired) + 1))
    
    print(f"Generating {total_reads} reads per file ({paired_percent}% paired)...")
    print(f"- Paired: {num_paired}")
    print(f"- Unpaired left: {len(left_unpaired_ids)}")
    print(f"- Unpaired right: {len(right_unpaired_ids)}")
    
    # Write left reads file
    with open(left_file_path, 'w') as left_file:
        # Write paired reads
        for read_id in paired_ids:
            sequence = generate_random_sequence(read_length)
            quality = generate_random_quality(read_length)
            write_fastq_record(left_file, generate_read_id(read_id, 1), sequence, quality)
        
        # Write unpaired reads
        for read_id in left_unpaired_ids:
            sequence = generate_random_sequence(read_length)
            quality = generate_random_quality(read_length)
            write_fastq_record(left_file, generate_read_id(read_id, 1), sequence, quality)
    
    # Write right reads file
    with open(right_file_path, 'w') as right_file:
        # Write paired reads
        for read_id in paired_ids:
            sequence = generate_random_sequence(read_length)
            quality = generate_random_quality(read_length)
            write_fastq_record(right_file, generate_read_id(read_id, 2), sequence, quality)
        
        # Write unpaired reads
        for read_id in right_unpaired_ids:
            sequence = generate_random_sequence(read_length)
            quality = generate_random_quality(read_length)
            write_fastq_record(right_file, generate_read_id(read_id, 2), sequence, quality)
    
    print(f"Created test datasets:")
    print(f"- Left file: {left_file_path}")
    print(f"- Right file: {right_file_path}")
    
    return left_file_path, right_file_path


def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Generate synthetic FASTQ datasets for testing get_pairs.py',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('-o', '--output-dir',
                        default='test_data',
                        help='Output directory for test datasets')
    
    parser.add_argument('-s', '--sizes',
                        nargs='+',
                        type=int,
                        default=[1000, 10000, 100000],
                        help='List of dataset sizes (reads per file) to generate')
    
    parser.add_argument('-p', '--paired-percents',
                        nargs='+',
                        type=float,
                        default=[50, 90, 100],
                        help='List of paired read percentages to generate')
    
    parser.add_argument('-l', '--read-length',
                        type=int,
                        default=150,
                        help='Length of each read')
    
    return parser.parse_args()


def main() -> int:
    """Main function."""
    args = parse_args()
    output_dir = Path(args.output_dir)
    
    # Create all combinations of sizes and paired percentages
    for size in args.sizes:
        for percent in args.paired_percents:
            dataset_dir = output_dir / f"size_{size}_paired_{int(percent)}"
            try:
                create_paired_fastq_files(
                    dataset_dir,
                    total_reads=size,
                    paired_percent=percent,
                    read_length=args.read_length
                )
            except Exception as e:
                print(f"Error creating dataset (size={size}, paired={percent}%): {e}")
                return 1
    
    return 0


if __name__ == '__main__':
    sys.exit(main()) 