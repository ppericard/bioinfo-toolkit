#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Get Pairs V2 Implementation

This module contains the improved implementation of get_pairs algorithm
using in-memory dictionaries for processing paired-end FASTQ files.
"""

from pathlib import Path
from typing import Dict, Set, Tuple

from bioinfotoolkit.utils.fastq_utils import open_file, extract_read_id
from bioinfotoolkit.scripts.fastq.get_pairs_common import (
    ensure_directory, count_reads, write_fastq_record, logger
)

class GetPairsV2:
    """Improved implementation of get_pairs using in-memory dictionaries."""
    
    @staticmethod
    def process(left_file: Path, right_file: Path, output_dir: Path, 
                compress: bool = False, verbose: bool = False) -> Dict[str, int]:
        """
        Process paired-end FASTQ files and separate them into pairs and singletons.
        
        Args:
            left_file: Path to the left (R1) FASTQ file
            right_file: Path to the right (R2) FASTQ file
            output_dir: Directory to write output files
            compress: Whether to compress output files
            verbose: Whether to print verbose output
            
        Returns:
            Dictionary with counts of paired and singleton reads
        """
        # Ensure output directory exists
        ensure_directory(output_dir)
        
        # Output file paths
        paired_1_path = output_dir / "paired_1.fastq"
        paired_2_path = output_dir / "paired_2.fastq"
        singleton_1_path = output_dir / "singleton_1.fastq"
        singleton_2_path = output_dir / "singleton_2.fastq"
        
        if compress:
            paired_1_path = paired_1_path.with_suffix(".fastq.gz")
            paired_2_path = paired_2_path.with_suffix(".fastq.gz")
            singleton_1_path = singleton_1_path.with_suffix(".fastq.gz")
            singleton_2_path = singleton_2_path.with_suffix(".fastq.gz")
        
        # Dictionaries to store reads
        left_reads = {}
        right_reads = {}
        
        # Process left file
        if verbose:
            logger.info(f"Processing left file: {left_file}")
        
        with open_file(left_file) as f:
            while True:
                header = f.readline().strip()
                if not header:
                    break
                
                sequence = f.readline().strip()
                plus_line = f.readline().strip()
                quality = f.readline().strip()
                
                read_id = extract_read_id(header)
                left_reads[read_id] = (header, sequence, plus_line, quality)
        
        # Process right file
        if verbose:
            logger.info(f"Processing right file: {right_file}")
        
        with open_file(right_file) as f:
            while True:
                header = f.readline().strip()
                if not header:
                    break
                
                sequence = f.readline().strip()
                plus_line = f.readline().strip()
                quality = f.readline().strip()
                
                read_id = extract_read_id(header)
                right_reads[read_id] = (header, sequence, plus_line, quality)
        
        # Find paired and singleton reads
        paired_ids = set(left_reads.keys()) & set(right_reads.keys())
        left_singleton_ids = set(left_reads.keys()) - paired_ids
        right_singleton_ids = set(right_reads.keys()) - paired_ids
        
        # Write paired reads
        if verbose:
            logger.info(f"Writing paired reads to {paired_1_path} and {paired_2_path}")
        
        with open_file(paired_1_path, 'w') as f1, open_file(paired_2_path, 'w') as f2:
            for read_id in paired_ids:
                header, sequence, plus_line, quality = left_reads[read_id]
                write_fastq_record(f1, header, sequence, plus_line, quality)
                
                header, sequence, plus_line, quality = right_reads[read_id]
                write_fastq_record(f2, header, sequence, plus_line, quality)
        
        # Write singleton reads
        if verbose:
            logger.info(f"Writing singleton reads to {singleton_1_path} and {singleton_2_path}")
        
        with open_file(singleton_1_path, 'w') as f:
            for read_id in left_singleton_ids:
                header, sequence, plus_line, quality = left_reads[read_id]
                write_fastq_record(f, header, sequence, plus_line, quality)
        
        with open_file(singleton_2_path, 'w') as f:
            for read_id in right_singleton_ids:
                header, sequence, plus_line, quality = right_reads[read_id]
                write_fastq_record(f, header, sequence, plus_line, quality)
        
        # Return counts
        counts = {
            'paired': len(paired_ids),
            'singleton_1': len(left_singleton_ids),
            'singleton_2': len(right_singleton_ids),
            'total_1': len(left_reads),
            'total_2': len(right_reads)
        }
        
        if verbose:
            logger.info(f"Paired reads: {counts['paired']}")
            logger.info(f"Singleton reads (left): {counts['singleton_1']}")
            logger.info(f"Singleton reads (right): {counts['singleton_2']}")
        
        return counts 