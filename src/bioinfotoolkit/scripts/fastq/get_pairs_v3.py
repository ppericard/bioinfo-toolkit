#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Get Pairs V3 Implementation

This module contains the memory-optimized implementation of get_pairs algorithm
using a disk-based approach for processing paired-end FASTQ files.
"""

import os
import shutil
import tempfile
from pathlib import Path
from typing import Dict, Set, Optional

from bioinfotoolkit.utils.fastq_utils import open_file, extract_read_id
from bioinfotoolkit.scripts.fastq.get_pairs_common import (
    ensure_directory, count_reads, write_fastq_record, logger
)

class GetPairsV3:
    """Memory-optimized implementation of get_pairs using a disk-based approach."""
    
    @staticmethod
    def process(left_file: Path, right_file: Path, output_dir: Path, 
                compress: bool = False, verbose: bool = False,
                chunk_size: int = 1000000, temp_dir: Optional[Path] = None) -> Dict[str, int]:
        """
        Process paired-end FASTQ files and separate them into pairs and singletons.
        
        This implementation uses a disk-based approach to minimize memory usage,
        making it suitable for processing very large FASTQ files.
        
        Args:
            left_file: Path to the left (R1) FASTQ file
            right_file: Path to the right (R2) FASTQ file
            output_dir: Directory to write output files
            compress: Whether to compress output files
            verbose: Whether to print verbose output
            chunk_size: Number of reads to process at once
            temp_dir: Directory for temporary files (default: system temp dir)
            
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
        
        # Create a temporary directory if not provided
        if temp_dir is None:
            temp_dir = Path(tempfile.mkdtemp())
            cleanup_temp = True
        else:
            ensure_directory(temp_dir)
            cleanup_temp = False
        
        try:
            # Step 1: Extract read IDs from both files
            if verbose:
                logger.info("Extracting read IDs...")
            
            left_ids_file = temp_dir / "left_ids.txt"
            right_ids_file = temp_dir / "right_ids.txt"
            
            # Extract and sort read IDs from left file
            with open(left_ids_file, 'w') as f:
                with open_file(left_file) as fastq:
                    for i, line in enumerate(fastq):
                        if i % 4 == 0:  # Header line
                            read_id = extract_read_id(line.strip())
                            f.write(f"{read_id}\n")
            
            # Extract and sort read IDs from right file
            with open(right_ids_file, 'w') as f:
                with open_file(right_file) as fastq:
                    for i, line in enumerate(fastq):
                        if i % 4 == 0:  # Header line
                            read_id = extract_read_id(line.strip())
                            f.write(f"{read_id}\n")
            
            # Step 2: Sort the ID files
            if verbose:
                logger.info("Sorting read IDs...")
            
            sorted_left_ids_file = temp_dir / "sorted_left_ids.txt"
            sorted_right_ids_file = temp_dir / "sorted_right_ids.txt"
            
            os.system(f"sort {left_ids_file} > {sorted_left_ids_file}")
            os.system(f"sort {right_ids_file} > {sorted_right_ids_file}")
            
            # Step 3: Find paired and singleton IDs
            if verbose:
                logger.info("Finding paired and singleton reads...")
            
            paired_ids_file = temp_dir / "paired_ids.txt"
            left_singleton_ids_file = temp_dir / "left_singleton_ids.txt"
            right_singleton_ids_file = temp_dir / "right_singleton_ids.txt"
            
            # Find paired IDs (intersection)
            os.system(f"comm -12 {sorted_left_ids_file} {sorted_right_ids_file} > {paired_ids_file}")
            
            # Find left singleton IDs (in left but not in right)
            os.system(f"comm -23 {sorted_left_ids_file} {sorted_right_ids_file} > {left_singleton_ids_file}")
            
            # Find right singleton IDs (in right but not in left)
            os.system(f"comm -13 {sorted_left_ids_file} {sorted_right_ids_file} > {right_singleton_ids_file}")
            
            # Step 4: Load IDs into memory
            with open(paired_ids_file) as f:
                paired_ids = set(line.strip() for line in f)
            
            with open(left_singleton_ids_file) as f:
                left_singleton_ids = set(line.strip() for line in f)
            
            with open(right_singleton_ids_file) as f:
                right_singleton_ids = set(line.strip() for line in f)
            
            # Step 5: Process and write output files
            if verbose:
                logger.info(f"Writing output files to {output_dir}...")
            
            # Process left file
            with open_file(paired_1_path, 'w') as paired_out, open_file(singleton_1_path, 'w') as singleton_out:
                with open_file(left_file) as f:
                    while True:
                        header = f.readline().strip()
                        if not header:
                            break
                        
                        sequence = f.readline().strip()
                        plus_line = f.readline().strip()
                        quality = f.readline().strip()
                        
                        read_id = extract_read_id(header)
                        
                        if read_id in paired_ids:
                            write_fastq_record(paired_out, header, sequence, plus_line, quality)
                        elif read_id in left_singleton_ids:
                            write_fastq_record(singleton_out, header, sequence, plus_line, quality)
            
            # Process right file
            with open_file(paired_2_path, 'w') as paired_out, open_file(singleton_2_path, 'w') as singleton_out:
                with open_file(right_file) as f:
                    while True:
                        header = f.readline().strip()
                        if not header:
                            break
                        
                        sequence = f.readline().strip()
                        plus_line = f.readline().strip()
                        quality = f.readline().strip()
                        
                        read_id = extract_read_id(header)
                        
                        if read_id in paired_ids:
                            write_fastq_record(paired_out, header, sequence, plus_line, quality)
                        elif read_id in right_singleton_ids:
                            write_fastq_record(singleton_out, header, sequence, plus_line, quality)
            
            # Return counts
            counts = {
                'paired': len(paired_ids),
                'singleton_1': len(left_singleton_ids),
                'singleton_2': len(right_singleton_ids),
                'total_1': len(paired_ids) + len(left_singleton_ids),
                'total_2': len(paired_ids) + len(right_singleton_ids)
            }
            
            if verbose:
                logger.info(f"Paired reads: {counts['paired']}")
                logger.info(f"Singleton reads (left): {counts['singleton_1']}")
                logger.info(f"Singleton reads (right): {counts['singleton_2']}")
            
            return counts
        
        finally:
            # Clean up temporary directory if we created it
            if cleanup_temp:
                shutil.rmtree(temp_dir) 