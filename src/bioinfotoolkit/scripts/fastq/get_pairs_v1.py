#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Get Pairs V1 Implementation

This module contains the original implementation of get_pairs algorithm (2012-2016)
for processing paired-end FASTQ files.
"""

import logging
from pathlib import Path
from typing import Dict, Set

from bioinfotoolkit.utils.fastq_utils import open_file, extract_read_id
from bioinfotoolkit.scripts.fastq.get_pairs_common import (
    ensure_directory, count_reads, build_read_id_set_from_file,
    get_output_paths, process_fastq_pairs, logger
)

class GetPairsV1:
    """
    Original implementation of get_pairs from 2012-2016.
    
    This is a simple implementation with minimal dependencies, using sets
    to find common read IDs between two files.
    """
    
    @staticmethod
    def process(left_file: Path, right_file: Path, output_dir: Path, 
                compress: bool = False, verbose: bool = False) -> Dict[str, int]:
        """
        Process paired-end FASTQ files to separate paired reads and singletons.
        
        Args:
            left_file: Path to the left reads FASTQ file
            right_file: Path to the right reads FASTQ file
            output_dir: Directory to write the output files
            compress: Whether to compress the output files
            verbose: Whether to print verbose progress information
            
        Returns:
            Dictionary with statistics about the processing
        """
        if verbose:
            logger.info("Using GetPairsV1 implementation (original, 2012-2016)")
            logger.info(f"Processing left file: {left_file}")
            logger.info(f"Processing right file: {right_file}")
            logger.info(f"Output directory: {output_dir}")
        
        # Ensure output directory exists
        ensure_directory(output_dir)
        
        # Get output file paths
        output_paths = get_output_paths(output_dir, compress)
        
        # Build read ID sets
        if verbose:
            logger.info("Building read ID sets...")
        
        left_ids = build_read_id_set_from_file(left_file, verbose)
        right_ids = build_read_id_set_from_file(right_file, verbose)
        
        # Find common read IDs
        common_ids = left_ids & right_ids
        
        if verbose:
            logger.info(f"Found {len(common_ids)} common reads")
            logger.info(f"Found {len(left_ids) - len(common_ids)} left singletons")
            logger.info(f"Found {len(right_ids) - len(common_ids)} right singletons")
        
        # Process and write output files
        if verbose:
            logger.info("Writing output files...")
        
        with open_file(output_paths['paired_1'], 'w') as paired_1_file, \
             open_file(output_paths['paired_2'], 'w') as paired_2_file, \
             open_file(output_paths['singleton_1'], 'w') as singleton_1_file, \
             open_file(output_paths['singleton_2'], 'w') as singleton_2_file:
            
            # Read and process left file
            left_singleton_count = 0
            paired_count = 0
            
            with open_file(left_file) as f:
                record = []
                for line_num, line in enumerate(f):
                    line = line.rstrip()
                    record.append(line)
                    
                    if (line_num + 1) % 4 == 0:  # End of FASTQ record
                        read_id = extract_read_id(record[0])
                        
                        if read_id in common_ids:
                            # Write to paired_1 file
                            paired_1_file.write('\n'.join(record) + '\n')
                            paired_count += 1
                        else:
                            # Write to singleton_1 file
                            singleton_1_file.write('\n'.join(record) + '\n')
                            left_singleton_count += 1
                        
                        record = []
            
            # Read and process right file
            right_singleton_count = 0
            
            with open_file(right_file) as f:
                record = []
                for line_num, line in enumerate(f):
                    line = line.rstrip()
                    record.append(line)
                    
                    if (line_num + 1) % 4 == 0:  # End of FASTQ record
                        read_id = extract_read_id(record[0])
                        
                        if read_id in common_ids:
                            # Write to paired_2 file
                            paired_2_file.write('\n'.join(record) + '\n')
                        else:
                            # Write to singleton_2 file
                            singleton_2_file.write('\n'.join(record) + '\n')
                            right_singleton_count += 1
                        
                        record = []
        
        # Return statistics
        return {
            'paired': paired_count,
            'left_singletons': left_singleton_count,
            'right_singletons': right_singleton_count
        } 