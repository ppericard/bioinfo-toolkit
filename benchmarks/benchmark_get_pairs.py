#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Benchmark Get Pairs

Description: Compare the performance of different implementations of get_pairs.py
             using synthetic test datasets.

This script benchmarks three implementations:
- v1: Original implementation from 2012-2016 (simple, minimal dependencies)
- v2: Improved implementation using in-memory dictionaries
- v3: Memory-optimized implementation using disk-based approach
"""

import os
import sys
import time
import argparse
import tempfile
import shutil
import json
import logging
import platform
import resource
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
import pandas as pd
import matplotlib.pyplot as plt

# Add the project root to the Python path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from tests.utils.test_data_generator import create_paired_fastq_files
from src.bioinfotoolkit.scripts.fastq.get_pairs_implementations import IMPLEMENTATIONS

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

def create_test_dataset(
    output_dir: Path,
    num_reads: int,
    paired_percent: float,
    read_length: int = 150
) -> Tuple[Path, Path]:
    """
    Create a test dataset with the specified parameters.
    
    Args:
        output_dir: Directory to write the dataset
        num_reads: Number of reads to generate
        paired_percent: Percentage of reads that should be paired (0-100)
        read_length: Length of each read
        
    Returns:
        Tuple of (left_file, right_file)
    """
    logger.info(f"Creating test dataset with {num_reads} reads ({paired_percent}% paired)")
    
    # Create the dataset directory
    dataset_dir = output_dir / f"reads_{num_reads}_{int(paired_percent)}"
    dataset_dir.mkdir(parents=True, exist_ok=True)
    
    # Generate the dataset
    left_file, right_file = create_paired_fastq_files(
        dataset_dir,
        total_reads=num_reads,
        paired_percent=paired_percent,
        read_length=read_length
    )
    
    logger.info(f"Created test dataset at {dataset_dir}")
    return left_file, right_file

def benchmark_implementation(
    implementation: str,
    left_file: Path,
    right_file: Path,
    output_dir: Path,
    compress: bool = False,
    verbose: bool = False,
    **kwargs
) -> Dict[str, Any]:
    """
    Benchmark a specific implementation of get_pairs.
    
    Args:
        implementation: Name of the implementation to benchmark
        left_file: Path to the left FASTQ file
        right_file: Path to the right FASTQ file
        output_dir: Directory to write output files
        compress: Whether to compress output files
        verbose: Whether to print verbose output
        **kwargs: Additional arguments to pass to the implementation
        
    Returns:
        Dictionary with benchmark results
    """
    logger.info(f"Benchmarking implementation: {implementation}")
    
    # Create output directory
    impl_output_dir = output_dir / implementation
    impl_output_dir.mkdir(parents=True, exist_ok=True)
    
    # Get the implementation class
    impl_class = IMPLEMENTATIONS[implementation]
    
    # Measure time and memory usage
    start_time = time.time()
    start_memory = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    
    # Run the implementation
    counts = impl_class.process(
        left_file,
        right_file,
        impl_output_dir,
        compress=compress,
        verbose=verbose,
        **kwargs
    )
    
    # Calculate elapsed time and memory usage
    elapsed_time = time.time() - start_time
    end_memory = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    memory_usage = end_memory - start_memory
    
    # Convert memory usage to MB (depends on platform)
    if platform.system() == 'Darwin':  # macOS
        memory_usage_mb = memory_usage / 1024 / 1024  # macOS reports in bytes
    else:  # Linux and others
        memory_usage_mb = memory_usage / 1024  # Linux reports in KB
    
    # Return benchmark results
    return {
        'implementation': implementation,
        'elapsed_time': elapsed_time,
        'memory_usage_mb': memory_usage_mb,
        'counts': counts
    }

def run_benchmarks(
    dataset_sizes: List[int],
    paired_percents: List[float],
    implementations: List[str],
    output_dir: Path,
    compress: bool = False,
    verbose: bool = False,
    read_length: int = 150,
    **kwargs
) -> List[Dict[str, Any]]:
    """
    Run benchmarks for all specified implementations and datasets.
    
    Args:
        dataset_sizes: List of dataset sizes to benchmark
        paired_percents: List of paired percentages to benchmark
        implementations: List of implementations to benchmark
        output_dir: Directory to write output files
        compress: Whether to compress output files
        verbose: Whether to print verbose output
        read_length: Length of each read
        **kwargs: Additional arguments to pass to the implementations
        
    Returns:
        List of benchmark results
    """
    results = []
    
    # Create datasets directory
    datasets_dir = output_dir / "datasets"
    datasets_dir.mkdir(parents=True, exist_ok=True)
    
    # Create results directory
    results_dir = output_dir / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    
    # Run benchmarks for each dataset size and paired percentage
    for size in dataset_sizes:
        for percent in paired_percents:
            # Create the test dataset
            left_file, right_file = create_test_dataset(
                datasets_dir,
                num_reads=size,
                paired_percent=percent,
                read_length=read_length
            )
            
            # Benchmark each implementation
            for implementation in implementations:
                try:
                    # Create a unique output directory for this benchmark
                    benchmark_output_dir = results_dir / f"size_{size}_paired_{int(percent)}"
                    benchmark_output_dir.mkdir(parents=True, exist_ok=True)
                    
                    # Run the benchmark
                    result = benchmark_implementation(
                        implementation,
                        left_file,
                        right_file,
                        benchmark_output_dir,
                        compress=compress,
                        verbose=verbose,
                        **kwargs
                    )
                    
                    # Add dataset information to the result
                    result['dataset_size'] = size
                    result['paired_percent'] = percent
                    
                    # Add the result to the list
                    results.append(result)
                    
                    # Log the result
                    logger.info(f"Implementation: {implementation}, "
                                f"Dataset: {size} reads ({percent}% paired), "
                                f"Time: {result['elapsed_time']:.2f}s, "
                                f"Memory: {result['memory_usage_mb']:.2f} MB")
                except Exception as e:
                    logger.error(f"Error benchmarking {implementation} on dataset "
                                f"size={size}, paired={percent}%: {e}")
    
    return results

def generate_report(results: List[Dict[str, Any]], output_dir: Path) -> None:
    """
    Generate a report from benchmark results.
    
    Args:
        results: List of benchmark results
        output_dir: Directory to write the report
    """
    # Convert results to a DataFrame
    df = pd.DataFrame(results)
    
    # Save the raw results as JSON
    results_file = output_dir / "benchmark_results.json"
    with open(results_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    # Save the DataFrame as CSV
    csv_file = output_dir / "benchmark_results.csv"
    df.to_csv(csv_file, index=False)
    
    # Generate plots
    generate_plots(df, output_dir)
    
    logger.info(f"Report generated at {output_dir}")

def generate_plots(df: pd.DataFrame, output_dir: Path) -> None:
    """
    Generate plots from benchmark results.
    
    Args:
        df: DataFrame with benchmark results
        output_dir: Directory to write the plots
    """
    # Create plots directory
    plots_dir = output_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    
    # Group by dataset size and implementation
    grouped = df.groupby(['dataset_size', 'implementation'])
    
    # Plot execution time by dataset size
    plt.figure(figsize=(10, 6))
    for name, group in grouped:
        size, impl = name
        plt.scatter(size, group['elapsed_time'].mean(), label=impl)
    
    plt.xlabel('Dataset Size (reads)')
    plt.ylabel('Execution Time (s)')
    plt.title('Execution Time by Dataset Size')
    plt.legend()
    plt.grid(True)
    plt.savefig(plots_dir / "execution_time_by_size.png")
    
    # Plot memory usage by dataset size
    plt.figure(figsize=(10, 6))
    for name, group in grouped:
        size, impl = name
        plt.scatter(size, group['memory_usage_mb'].mean(), label=impl)
    
    plt.xlabel('Dataset Size (reads)')
    plt.ylabel('Memory Usage (MB)')
    plt.title('Memory Usage by Dataset Size')
    plt.legend()
    plt.grid(True)
    plt.savefig(plots_dir / "memory_usage_by_size.png")
    
    # Plot execution time by implementation
    plt.figure(figsize=(12, 6))
    implementations = df['implementation'].unique()
    dataset_sizes = df['dataset_size'].unique()
    
    x = np.arange(len(implementations))
    width = 0.8 / len(dataset_sizes)
    
    for i, size in enumerate(dataset_sizes):
        times = [df[(df['implementation'] == impl) & (df['dataset_size'] == size)]['elapsed_time'].mean() 
                for impl in implementations]
        plt.bar(x + i * width, times, width, label=f'{size} reads')
    
    plt.xlabel('Implementation')
    plt.ylabel('Execution Time (s)')
    plt.title('Execution Time by Implementation')
    plt.xticks(x + width * (len(dataset_sizes) - 1) / 2, implementations)
    plt.legend()
    plt.grid(True)
    plt.savefig(plots_dir / "execution_time_by_implementation.png")
    
    # Plot memory usage by implementation
    plt.figure(figsize=(12, 6))
    
    for i, size in enumerate(dataset_sizes):
        memory = [df[(df['implementation'] == impl) & (df['dataset_size'] == size)]['memory_usage_mb'].mean() 
                 for impl in implementations]
        plt.bar(x + i * width, memory, width, label=f'{size} reads')
    
    plt.xlabel('Implementation')
    plt.ylabel('Memory Usage (MB)')
    plt.title('Memory Usage by Implementation')
    plt.xticks(x + width * (len(dataset_sizes) - 1) / 2, implementations)
    plt.legend()
    plt.grid(True)
    plt.savefig(plots_dir / "memory_usage_by_implementation.png")

def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Benchmark different implementations of get_pairs",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        '-o', '--output-dir',
        default='benchmark_results',
        help='Output directory for benchmark results'
    )
    
    parser.add_argument(
        '-s', '--sizes',
        nargs='+',
        type=int,
        default=[1000, 10000, 100000],
        help='List of dataset sizes to benchmark'
    )
    
    parser.add_argument(
        '-p', '--paired-percents',
        nargs='+',
        type=float,
        default=[50, 90, 100],
        help='List of paired percentages to benchmark'
    )
    
    parser.add_argument(
        '-i', '--implementations',
        nargs='+',
        choices=list(IMPLEMENTATIONS.keys()),
        default=list(IMPLEMENTATIONS.keys()),
        help='List of implementations to benchmark'
    )
    
    parser.add_argument(
        '-l', '--read-length',
        type=int,
        default=150,
        help='Length of each read in the test datasets'
    )
    
    parser.add_argument(
        '-z', '--compress',
        action='store_true',
        help='Compress output files with gzip'
    )
    
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Print verbose output'
    )
    
    # V3-specific options
    parser.add_argument(
        '--chunk-size',
        type=int,
        default=1000000,
        help='Number of reads to process at once (V3 only)'
    )
    
    parser.add_argument(
        '--temp-dir',
        help='Directory for temporary files (V3 only)'
    )
    
    return parser.parse_args()

def main() -> int:
    """Main function."""
    # Parse command line arguments
    args = parse_args()
    
    try:
        # Create output directory
        output_dir = Path(args.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Prepare kwargs for the implementations
        kwargs = {}
        if 'v3' in args.implementations:
            kwargs['chunk_size'] = args.chunk_size
            if args.temp_dir:
                kwargs['temp_dir'] = Path(args.temp_dir)
        
        # Run benchmarks
        results = run_benchmarks(
            dataset_sizes=args.sizes,
            paired_percents=args.paired_percents,
            implementations=args.implementations,
            output_dir=output_dir,
            compress=args.compress,
            verbose=args.verbose,
            read_length=args.read_length,
            **kwargs
        )
        
        # Generate report
        generate_report(results, output_dir)
        
        return 0
    except Exception as e:
        logger.error(f"Error: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    # Import numpy here to avoid importing it in the module scope
    import numpy as np
    sys.exit(main()) 