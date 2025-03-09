#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Benchmark Get Pairs

Description: Compare the performance of the original and new implementations of get_pairs.py
             using synthetic test datasets.
"""

import os
import sys
import time
import argparse
import subprocess
import tempfile
import shutil
from pathlib import Path
import resource
import json
import pandas as pd
import matplotlib.pyplot as plt
from typing import Dict, List, Tuple, Optional

# Path to implementations
ORIGINAL_SCRIPT = 'bin/get_pairs.py.original'
NEW_SCRIPT = 'bin/get_pairs.py'


def copy_original_implementation():
    """Make a copy of the current get_pairs.py as the original implementation."""
    if not os.path.exists(ORIGINAL_SCRIPT):
        shutil.copy(NEW_SCRIPT, ORIGINAL_SCRIPT)
        print(f"Copied current implementation to {ORIGINAL_SCRIPT}")


def create_benchmark_datasets(sizes: List[int], paired_percents: List[float], output_dir: str) -> None:
    """Create the test datasets using the create_test_datasets.py script."""
    cmd = [
        'python', 'bin/create_test_datasets.py',
        '-o', output_dir,
        '-s', *[str(s) for s in sizes],
        '-p', *[str(p) for p in paired_percents]
    ]
    
    print(f"Creating benchmark datasets: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)


def run_benchmark(
    test_data_dir: str,
    results_file: str,
    sizes: List[int],
    paired_percents: List[float],
    repeats: int = 3
) -> None:
    """
    Run the benchmark comparing original and new implementations.
    
    Args:
        test_data_dir: Directory containing test datasets
        results_file: Path to write results as JSON
        sizes: List of dataset sizes to test
        paired_percents: List of paired percentages to test
        repeats: Number of times to repeat each test
    """
    results = []
    test_data_dir = Path(test_data_dir)
    
    # Loop through all test combinations
    for size in sizes:
        for percent in paired_percents:
            dataset_dir = test_data_dir / f"size_{size}_paired_{int(percent)}"
            left_path = dataset_dir / f"reads_1_{size}_{int(percent)}.fastq"
            right_path = dataset_dir / f"reads_2_{size}_{int(percent)}.fastq"
            
            if not left_path.exists() or not right_path.exists():
                print(f"Warning: Test dataset not found at {dataset_dir}")
                continue
            
            print(f"\nBenchmarking: size={size}, paired={percent}%")
            
            # Repeat the test multiple times for both implementations
            for implementation in ['original', 'new']:
                script = ORIGINAL_SCRIPT if implementation == 'original' else NEW_SCRIPT
                
                if not os.path.exists(script):
                    print(f"Warning: Implementation {script} not found, skipping.")
                    continue
                
                for run in range(1, repeats + 1):
                    print(f"  Running {implementation} implementation (run {run}/{repeats})...")
                    
                    # Create a temporary output directory
                    with tempfile.TemporaryDirectory() as temp_dir:
                        cmd = [
                            'python', script,
                            '-l', str(left_path),
                            '-r', str(right_path),
                            '-o', temp_dir
                        ]
                        
                        # Measure time and memory usage
                        start_time = time.time()
                        
                        # Use resource module to get max memory usage
                        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                        process.wait()
                        
                        elapsed_time = time.time() - start_time
                        
                        # Get memory usage by reading output files (estimate based on file sizes)
                        output_size = sum(os.path.getsize(os.path.join(temp_dir, f)) for f in os.listdir(temp_dir))
                        
                        # Alternative: use /usr/bin/time -v to get peak memory usage
                        try:
                            cmd_time = ['/usr/bin/time', '-v', 'python', script, 
                                       '-l', str(left_path), 
                                       '-r', str(right_path), 
                                       '-o', temp_dir]
                            time_output = subprocess.check_output(cmd_time, stderr=subprocess.STDOUT, text=True)
                            
                            # Parse time output to get max memory
                            max_memory = 0
                            for line in time_output.splitlines():
                                if 'Maximum resident set size' in line:
                                    max_memory = int(line.split(':')[1].strip())
                                    break
                        except (subprocess.SubprocessError, FileNotFoundError):
                            # If /usr/bin/time fails, use a fallback estimation
                            max_memory = output_size * 2  # Simple estimation
                        
                        # Store results
                        result = {
                            'implementation': implementation,
                            'size': size,
                            'paired_percent': percent,
                            'run': run,
                            'time_seconds': elapsed_time,
                            'max_memory_kb': max_memory / 1024,  # Convert to KB
                            'output_size_bytes': output_size
                        }
                        
                        results.append(result)
                        print(f"    Time: {elapsed_time:.2f}s, Memory: {max_memory/1024/1024:.2f}MB")
    
    # Save results
    with open(results_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\nResults saved to {results_file}")


def analyze_results(results_file: str, output_dir: str) -> None:
    """
    Analyze benchmark results and create visualizations.
    
    Args:
        results_file: Path to JSON results file
        output_dir: Directory to save plots and reports
    """
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Load results
    with open(results_file, 'r') as f:
        results = json.load(f)
    
    # Convert to pandas DataFrame
    df = pd.DataFrame(results)
    
    # Group by implementation, size, and paired_percent, then calculate mean and std
    summary = df.groupby(['implementation', 'size', 'paired_percent']).agg({
        'time_seconds': ['mean', 'std'],
        'max_memory_kb': ['mean', 'std']
    }).reset_index()
    
    # Create plots
    # 1. Time vs Size by Implementation
    plt.figure(figsize=(10, 6))
    for impl in df['implementation'].unique():
        data = summary[summary['implementation'] == impl]
        plt.errorbar(
            data['size'], 
            data['time_seconds']['mean'],
            yerr=data['time_seconds']['std'],
            label=f"{impl.capitalize()} Implementation",
            marker='o'
        )
    
    plt.xlabel('Dataset Size (reads per file)')
    plt.ylabel('Execution Time (seconds)')
    plt.title('Execution Time vs Dataset Size')
    plt.legend()
    plt.grid(True)
    plt.xscale('log')
    plt.savefig(os.path.join(output_dir, 'time_vs_size.png'))
    
    # 2. Memory vs Size by Implementation
    plt.figure(figsize=(10, 6))
    for impl in df['implementation'].unique():
        data = summary[summary['implementation'] == impl]
        plt.errorbar(
            data['size'], 
            data['max_memory_kb']['mean'] / 1024,  # Convert to MB
            yerr=data['max_memory_kb']['std'] / 1024,
            label=f"{impl.capitalize()} Implementation",
            marker='o'
        )
    
    plt.xlabel('Dataset Size (reads per file)')
    plt.ylabel('Memory Usage (MB)')
    plt.title('Memory Usage vs Dataset Size')
    plt.legend()
    plt.grid(True)
    plt.xscale('log')
    plt.savefig(os.path.join(output_dir, 'memory_vs_size.png'))
    
    # 3. Time Ratio (Original/New) vs Size
    plt.figure(figsize=(10, 6))
    time_ratio = []
    size_list = []
    
    for size in df['size'].unique():
        for percent in df['paired_percent'].unique():
            orig = summary[(summary['implementation'] == 'original') & 
                          (summary['size'] == size) & 
                          (summary['paired_percent'] == percent)]
            new = summary[(summary['implementation'] == 'new') & 
                         (summary['size'] == size) & 
                         (summary['paired_percent'] == percent)]
            
            if not orig.empty and not new.empty:
                ratio = orig['time_seconds']['mean'].values[0] / new['time_seconds']['mean'].values[0]
                time_ratio.append(ratio)
                size_list.append(size)
    
    plt.scatter(size_list, time_ratio, s=80)
    plt.axhline(y=1.0, color='r', linestyle='--', label='Equal Performance')
    plt.xlabel('Dataset Size (reads per file)')
    plt.ylabel('Time Ratio (Original/New)')
    plt.title('Performance Improvement Ratio vs Dataset Size')
    plt.grid(True)
    plt.xscale('log')
    plt.legend()
    plt.savefig(os.path.join(output_dir, 'time_ratio.png'))
    
    # Save summary to CSV
    summary.to_csv(os.path.join(output_dir, 'summary.csv'))
    
    # Create a text report
    with open(os.path.join(output_dir, 'report.txt'), 'w') as f:
        f.write("Get Pairs Performance Benchmark Report\n")
        f.write("====================================\n\n")
        
        f.write("Summary Statistics:\n")
        f.write("-----------------\n")
        
        # Average improvement
        avg_time_ratio = sum(time_ratio) / len(time_ratio) if time_ratio else 0
        f.write(f"Average Speed Improvement: {avg_time_ratio:.2f}x faster\n\n")
        
        # Detailed results
        f.write("Detailed Results:\n")
        f.write("----------------\n")
        
        for size in sorted(df['size'].unique()):
            f.write(f"\nDataset Size: {size} reads per file\n")
            
            for percent in sorted(df['paired_percent'].unique()):
                f.write(f"  Paired Percentage: {percent}%\n")
                
                orig = df[(df['implementation'] == 'original') & 
                         (df['size'] == size) & 
                         (df['paired_percent'] == percent)]
                
                new = df[(df['implementation'] == 'new') & 
                        (df['size'] == size) & 
                        (df['paired_percent'] == percent)]
                
                if not orig.empty and not new.empty:
                    orig_time = orig['time_seconds'].mean()
                    new_time = new['time_seconds'].mean()
                    
                    orig_mem = orig['max_memory_kb'].mean() / 1024  # Convert to MB
                    new_mem = new['max_memory_kb'].mean() / 1024
                    
                    time_ratio = orig_time / new_time if new_time > 0 else float('inf')
                    mem_ratio = orig_mem / new_mem if new_mem > 0 else float('inf')
                    
                    f.write(f"    Original: {orig_time:.2f}s, {orig_mem:.2f}MB\n")
                    f.write(f"    New:      {new_time:.2f}s, {new_mem:.2f}MB\n")
                    f.write(f"    Improvement: {time_ratio:.2f}x faster, {mem_ratio:.2f}x memory efficient\n")
    
    print(f"Analysis and visualizations saved to {output_dir}")


def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Benchmark get_pairs.py implementations',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('-o', '--output-dir',
                        default='benchmark_results',
                        help='Output directory for benchmark results')
    
    parser.add_argument('-d', '--data-dir',
                        default='test_data',
                        help='Directory for test datasets')
    
    parser.add_argument('-s', '--sizes',
                        nargs='+',
                        type=int,
                        default=[1000, 10000, 100000],
                        help='List of dataset sizes to benchmark')
    
    parser.add_argument('-p', '--paired-percents',
                        nargs='+',
                        type=float,
                        default=[50, 90, 100],
                        help='List of paired read percentages to benchmark')
    
    parser.add_argument('-r', '--repeats',
                        type=int,
                        default=3,
                        help='Number of times to repeat each test')
    
    parser.add_argument('-c', '--create-datasets',
                        action='store_true',
                        help='Create test datasets before benchmarking')
    
    return parser.parse_args()


def main() -> int:
    """Main function."""
    args = parse_args()
    output_dir = Path(args.output_dir)
    data_dir = Path(args.data_dir)
    results_file = output_dir / 'results.json'
    
    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Check if both implementations exist
    if not os.path.exists(NEW_SCRIPT):
        print(f"Error: New implementation '{NEW_SCRIPT}' not found")
        return 1
    
    # Create a copy of the current implementation as the original
    copy_original_implementation()
    
    # Create test datasets if requested
    if args.create_datasets:
        create_benchmark_datasets(args.sizes, args.paired_percents, str(data_dir))
    
    # Run benchmarks
    try:
        run_benchmark(
            test_data_dir=str(data_dir),
            results_file=str(results_file),
            sizes=args.sizes,
            paired_percents=args.paired_percents,
            repeats=args.repeats
        )
    except Exception as e:
        print(f"Error running benchmarks: {e}")
        return 1
    
    # Analyze results
    try:
        analyze_results(
            results_file=str(results_file),
            output_dir=str(output_dir / 'analysis')
        )
    except Exception as e:
        print(f"Error analyzing results: {e}")
        return 1
    
    return 0


if __name__ == '__main__':
    sys.exit(main()) 