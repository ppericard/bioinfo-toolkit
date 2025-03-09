#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Create Test Datasets

Command-line utility to generate synthetic datasets for testing.
"""

import os
import sys
import argparse
from pathlib import Path

from tests.utils.test_data_generator import (
    create_paired_fastq_files,
    create_fasta_dataset
)

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate synthetic datasets for testing",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        '-o', '--output-dir',
        default='tests/data/generated',
        help='Output directory for test datasets'
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Dataset type to generate')
    subparsers.required = True
    
    # FASTQ paired-end datasets
    fastq_parser = subparsers.add_parser('fastq', help='Generate paired-end FASTQ datasets')
    fastq_parser.add_argument(
        '-n', '--num-reads',
        type=int,
        default=1000,
        help='Number of reads to generate'
    )
    fastq_parser.add_argument(
        '-p', '--paired-percent',
        type=float,
        default=90.0,
        help='Percentage of reads that should be paired (0-100)'
    )
    fastq_parser.add_argument(
        '-l', '--read-length',
        type=int,
        default=150,
        help='Length of each read'
    )
    fastq_parser.add_argument(
        '-d', '--dataset-name',
        default='test_dataset',
        help='Name of the dataset (will be used as a subdirectory)'
    )
    
    # FASTA datasets
    fasta_parser = subparsers.add_parser('fasta', help='Generate FASTA datasets')
    fasta_parser.add_argument(
        '-n', '--num-sequences',
        type=int,
        default=100,
        help='Number of sequences to generate'
    )
    fasta_parser.add_argument(
        '--min-length',
        type=int,
        default=100,
        help='Minimum sequence length'
    )
    fasta_parser.add_argument(
        '--max-length',
        type=int,
        default=1000,
        help='Maximum sequence length'
    )
    fasta_parser.add_argument(
        '-d', '--dataset-name',
        default='test_dataset',
        help='Name of the dataset (will be used for filename)'
    )
    
    return parser.parse_args()

def main():
    """Main function."""
    args = parse_args()
    output_dir = Path(args.output_dir)
    
    if args.command == 'fastq':
        dataset_dir = output_dir / args.dataset_name
        try:
            r1_file, r2_file = create_paired_fastq_files(
                dataset_dir,
                total_reads=args.num_reads,
                paired_percent=args.paired_percent,
                read_length=args.read_length
            )
            print(f"Created paired FASTQ dataset:")
            print(f"- R1 file: {r1_file}")
            print(f"- R2 file: {r2_file}")
            print(f"- Total reads: {args.num_reads}")
            print(f"- Paired reads: {int(args.num_reads * args.paired_percent / 100)}")
        except Exception as e:
            print(f"Error creating FASTQ dataset: {e}")
            return 1
    
    elif args.command == 'fasta':
        try:
            output_path = output_dir / f"{args.dataset_name}.fasta"
            fasta_file = create_fasta_dataset(
                output_path,
                num_sequences=args.num_sequences,
                min_length=args.min_length,
                max_length=args.max_length
            )
            print(f"Created FASTA dataset:")
            print(f"- File: {fasta_file}")
            print(f"- Sequences: {args.num_sequences}")
            print(f"- Length range: {args.min_length}-{args.max_length}")
        except Exception as e:
            print(f"Error creating FASTA dataset: {e}")
            return 1
    
    return 0

if __name__ == "__main__":
    sys.exit(main()) 