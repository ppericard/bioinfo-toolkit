#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Test Get Pairs

Description: Validate and benchmark the get_pairs.py implementations
"""

import os
import sys
import time
import argparse
import tempfile
import shutil
import logging
from pathlib import Path
import subprocess
import json
import platform
from typing import Dict, Any, Tuple, Optional, List, Union

try:
    from bioinfo_logger import get_logger
    from memory_tracker import measure_script_memory, format_bytes, MemoryTracker
except ImportError:
    sys.path.append(str(Path(__file__).parent))
    from bioinfo_logger import get_logger
    from memory_tracker import measure_script_memory, format_bytes, MemoryTracker

# Constants
ORIGINAL_SCRIPT = 'bin/get_pairs.py.original'
V2_SCRIPT = 'bin/get_pairs.py'
V3_SCRIPT = 'bin/get_pairs_v3.py'
DEFAULT_TEST_SIZES = [1000, 10000, 100000]
DEFAULT_PAIRED_PERCENTS = [50, 90, 100]


def create_test_datasets(test_dir: str, sizes: List[int], paired_percents: List[float]) -> dict:
    """
    Create test datasets using create_test_datasets.py.
    
    Args:
        test_dir: Directory to store test datasets
        sizes: List of dataset sizes to create
        paired_percents: List of paired percentages to create
        
    Returns:
        Dictionary mapping test cases to file paths
    """
    logger = get_logger("test_get_pairs")
    logger.info("Creating test datasets")
    
    # Ensure test directory exists
    test_dir_path = Path(test_dir)
    test_dir_path.mkdir(parents=True, exist_ok=True)
    
    # Run create_test_datasets.py script
    sizes_str = " ".join(str(s) for s in sizes)
    percents_str = " ".join(str(p) for p in paired_percents)
    
    cmd = [
        sys.executable, 'bin/create_test_datasets.py',
        '-o', test_dir,
        '-s', *sizes_str.split(),
        '-p', *percents_str.split()
    ]
    
    logger.info(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, check=True)
    
    # Build a dictionary of test case paths
    test_cases = {}
    for size in sizes:
        for percent in paired_percents:
            case_dir = test_dir_path / f"size_{size}_paired_{int(percent)}"
            left_file = case_dir / f"reads_1_{size}_{int(percent)}.fastq"
            right_file = case_dir / f"reads_2_{size}_{int(percent)}.fastq"
            
            if left_file.exists() and right_file.exists():
                test_cases[f"{size}_{percent}"] = {
                    "left": str(left_file),
                    "right": str(right_file),
                    "size": size,
                    "paired_percent": percent
                }
    
    logger.info(f"Created {len(test_cases)} test cases")
    return test_cases


def get_directory_size(path: Union[str, Path]) -> int:
    """
    Get the total size of a directory in bytes.
    
    Args:
        path: Directory path
        
    Returns:
        Size in bytes
    """
    total_size = 0
    path = Path(path)
    
    for entry in path.glob('**/*'):
        if entry.is_file():
            total_size += entry.stat().st_size
            
    return total_size


def run_get_pairs(
    script_path: str, 
    left_file: str, 
    right_file: str, 
    output_dir: str, 
    verbose: bool = False,
    keep_temp: bool = False
) -> Tuple[float, Dict[str, Any], Dict[str, Any], int]:
    """
    Run get_pairs script, measuring execution time, memory usage, and disk usage.
    
    Args:
        script_path: Path to get_pairs.py script
        left_file: Path to left reads file
        right_file: Path to right reads file
        output_dir: Output directory
        verbose: Enable verbose output
        keep_temp: Keep temporary files
        
    Returns:
        Tuple of (execution time in seconds, memory statistics, output file stats, disk usage)
    """
    # Create a temp directory for measuring disk usage
    temp_dir = None
    if "v3.py" in script_path:
        temp_dir = tempfile.mkdtemp(prefix="get_pairs_test_")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Build command
    cmd = [
        sys.executable, script_path,
        '-l', left_file,
        '-r', right_file,
        '-o', output_dir
    ]
    
    if verbose:
        cmd.append('-v')
        
    if "v3.py" in script_path and temp_dir:
        cmd.extend(['--temp-dir', temp_dir])
        if keep_temp:
            cmd.append('--keep-temp')
    
    # Measure time and memory usage
    start_time = time.time()
    output, memory_stats = measure_script_memory(cmd)
    end_time = time.time()
    
    # Calculate execution time
    execution_time = end_time - start_time
    
    # Get disk usage
    disk_usage = 0
    if temp_dir and os.path.exists(temp_dir):
        disk_usage = get_directory_size(temp_dir)
        if not keep_temp:
            shutil.rmtree(temp_dir, ignore_errors=True)
    
    # Get statistics from output files
    output_stats = {
        "left_paired": 0,
        "left_unpaired": 0,
        "right_paired": 0,
        "right_unpaired": 0
    }
    
    # Count records in output files
    output_dir_path = Path(output_dir)
    for direction in ["left", "right"]:
        for pair_type in ["paired", "unpaired"]:
            # Look for output files with various extensions
            found = False
            for ext in ["fastq", "fq", "fastq.gz", "fq.gz"]:
                pattern = f"*.{pair_type}.{ext}"
                matching_files = list(output_dir_path.glob(pattern))
                
                if matching_files:
                    file_path = matching_files[0]
                    # Count records (every 4 lines is a record in FASTQ)
                    line_count = 0
                    with open(file_path, 'r') as f:
                        for line in f:
                            line_count += 1
                    
                    output_stats[f"{direction}_{pair_type}"] = line_count // 4
                    found = True
                    break
            
            if not found:
                print(f"Warning: No {pair_type} file found for {direction} reads")
    
    return execution_time, memory_stats, output_stats, disk_usage


def test_get_pairs(
    test_cases: Dict[str, Dict[str, Any]], 
    repeats: int = 3, 
    save_path: Optional[str] = None,
    selected_implementations: Optional[List[str]] = None
) -> Dict[str, Any]:
    """
    Test all three implementations of get_pairs.py and compare performance.
    
    Args:
        test_cases: Dictionary mapping test case names to file paths
        repeats: Number of times to repeat each test
        save_path: Optional path to save results as JSON
        selected_implementations: Optional list of implementation names to test
        
    Returns:
        Dictionary with benchmark results
    """
    logger = get_logger("test_get_pairs")
    results = {}
    
    # Check which implementations are available
    all_implementations = [
        ("original", ORIGINAL_SCRIPT),
        ("v2", V2_SCRIPT),
        ("v3", V3_SCRIPT)
    ]
    
    # Filter implementations that exist and are selected
    implementations = []
    for impl_name, script_path in all_implementations:
        if os.path.exists(script_path):
            if selected_implementations is None or impl_name in selected_implementations:
                implementations.append((impl_name, script_path))
    
    if not implementations:
        logger.error("No implementations found")
        return {}
    
    logger.info(f"Testing implementations: {', '.join(impl[0] for impl in implementations)}")
    
    # Run tests for each test case
    for case_name, case_info in test_cases.items():
        logger.info(f"Testing case: {case_name}")
        results[case_name] = {
            impl_name: {
                "times": [], 
                "memory": [], 
                "disk_usage": [],
                "stats": {}
            }
            for impl_name, _ in implementations
        }
        results[case_name]["case_info"] = case_info
        
        # Test all implementations
        for impl_name, script_path in implementations:
            logger.info(f"  Running {impl_name} implementation...")
            
            # Run multiple times for reliable results
            for run in range(1, repeats + 1):
                logger.info(f"    Run {run}/{repeats}")
                
                # Create a temporary output directory for this run
                with tempfile.TemporaryDirectory() as temp_dir:
                    try:
                        # Run the script and measure time, memory, and disk usage
                        execution_time, memory_stats, output_stats, disk_usage = run_get_pairs(
                            script_path=script_path,
                            left_file=case_info["left"],
                            right_file=case_info["right"],
                            output_dir=temp_dir,
                            verbose=False,
                            keep_temp=False
                        )
                        
                        # Store results
                        results[case_name][impl_name]["times"].append(execution_time)
                        results[case_name][impl_name]["memory"].append(memory_stats.get('peak_bytes', 0))
                        results[case_name][impl_name]["disk_usage"].append(disk_usage)
                        
                        # Only store output stats from the last run
                        results[case_name][impl_name]["stats"] = output_stats
                        
                        # Format memory and disk usage for display
                        peak_memory = memory_stats.get('peak_bytes', 0)
                        formatted_memory = format_bytes(peak_memory)
                        formatted_disk = format_bytes(disk_usage)
                        
                        logger.info(f"      Time: {execution_time:.4f}s, Memory: {formatted_memory}, Disk: {formatted_disk}")
                    except Exception as e:
                        logger.error(f"Error running {impl_name} implementation: {e}")
                        if hasattr(e, '__traceback__'):
                            import traceback
                            logger.error(traceback.format_exc())
    
    # Save results to JSON if requested
    if save_path:
        # Convert results to serializable format (remove non-JSON serializable types)
        serializable_results = {
            case: {
                impl: {
                    "times": results[case][impl]["times"],
                    "memory": results[case][impl]["memory"],
                    "disk_usage": results[case][impl]["disk_usage"],
                    "stats": results[case][impl]["stats"],
                }
                for impl in results[case] if impl != "case_info"
            }
            for case in results
        }
        
        # Add metadata
        serializable_results["metadata"] = {
            "platform": platform.system(),
            "python_version": platform.python_version(),
            "date": time.strftime("%Y-%m-%d %H:%M:%S"),
            "repeats": repeats
        }
        
        # Save to file
        with open(save_path, 'w') as f:
            json.dump(serializable_results, f, indent=2)
        
        logger.info(f"Results saved to {save_path}")
    
    return results


def analyze_results(results: Dict[str, Any]) -> None:
    """
    Analyze benchmark results and print summary.
    
    Args:
        results: Dictionary with benchmark results
    """
    logger = get_logger("test_get_pairs")
    logger.info("\nBenchmark Results:")
    logger.info("=================")
    
    # Get all implementations
    implementations = []
    for case_name, case_results in results.items():
        if case_name == "metadata":
            continue
        implementations = [impl for impl in case_results if impl != "case_info"]
        if implementations:
            break
    
    if not implementations:
        logger.error("No implementation data found")
        return
    
    # Calculate averages and improvements
    for case_name, case_results in results.items():
        if case_name == "metadata":
            continue
            
        logger.info(f"\nCase: {case_name}")
        
        # Skip if no data for any implementation
        if not any(case_results.get(impl, {}).get("times") for impl in implementations):
            logger.warning("  No data available")
            continue
        
        # Use the first implementation as a reference for comparisons
        reference_impl = implementations[0]
        
        # Calculate timing statistics for all implementations
        logger.info("  Execution Time:")
        impl_time_stats = {}
        
        for impl in implementations:
            times = case_results.get(impl, {}).get("times", [])
            if not times:
                continue
                
            avg_time = sum(times) / len(times)
            std_time = (sum((t - avg_time) ** 2 for t in times) / len(times)) ** 0.5
            
            impl_time_stats[impl] = {
                "avg": avg_time,
                "std": std_time
            }
            
            logger.info(f"    {impl.ljust(10)}: {avg_time:.4f}s ± {std_time:.4f}s")
        
        # Calculate memory usage statistics
        logger.info("  Memory Usage:")
        impl_memory_stats = {}
        
        for impl in implementations:
            memory_values = case_results.get(impl, {}).get("memory", [])
            if not memory_values:
                continue
                
            avg_memory = sum(memory_values) / len(memory_values)
            std_memory = (sum((m - avg_memory) ** 2 for m in memory_values) / len(memory_values)) ** 0.5
            
            impl_memory_stats[impl] = {
                "avg": avg_memory,
                "std": std_memory
            }
            
            logger.info(f"    {impl.ljust(10)}: {format_bytes(int(avg_memory))} ± {format_bytes(int(std_memory))}")
        
        # Calculate disk usage statistics for V3
        if "v3" in implementations:
            logger.info("  Disk Usage:")
            disk_values = case_results.get("v3", {}).get("disk_usage", [])
            if disk_values:
                avg_disk = sum(disk_values) / len(disk_values)
                std_disk = (sum((d - avg_disk) ** 2 for d in disk_values) / len(disk_values)) ** 0.5
                
                logger.info(f"    {'v3'.ljust(10)}: {format_bytes(int(avg_disk))} ± {format_bytes(int(std_disk))}")
        
        # Compare implementations to the reference
        logger.info("  Performance Comparison:")
        
        for impl in implementations[1:]:
            if impl not in impl_time_stats or reference_impl not in impl_time_stats:
                continue
                
            # Time comparison
            time_ratio = impl_time_stats[reference_impl]["avg"] / impl_time_stats[impl]["avg"]
            if time_ratio > 1:
                logger.info(f"    {impl} is {time_ratio:.2f}x faster than {reference_impl}")
            else:
                logger.info(f"    {impl} is {1/time_ratio:.2f}x slower than {reference_impl}")
            
            # Memory comparison
            if impl in impl_memory_stats and reference_impl in impl_memory_stats:
                memory_ratio = impl_memory_stats[reference_impl]["avg"] / impl_memory_stats[impl]["avg"]
                if memory_ratio > 1:
                    logger.info(f"    {impl} uses {memory_ratio:.2f}x less memory than {reference_impl}")
                else:
                    logger.info(f"    {impl} uses {1/memory_ratio:.2f}x more memory than {reference_impl}")
        
        # Validate output stats
        logger.info("  Output Validation:")
        
        reference_stats = case_results.get(reference_impl, {}).get("stats", {})
        validation_passed = True
        
        for impl in implementations[1:]:
            impl_stats = case_results.get(impl, {}).get("stats", {})
            if not impl_stats or not reference_stats:
                continue
                
            mismatches = []
            for key in reference_stats:
                if key in impl_stats and reference_stats[key] != impl_stats[key]:
                    mismatches.append(f"{key}: {reference_stats[key]} vs {impl_stats[key]}")
            
            if mismatches:
                logger.warning(f"    {impl} vs {reference_impl}: Differences detected: {', '.join(mismatches)}")
                validation_passed = False
        
        if validation_passed:
            logger.info("    PASSED - All implementations produced identical results")


def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Test and benchmark get_pairs.py implementations',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('-t', '--test-dir',
                        default='test_data',
                        help='Directory for test datasets')
    
    parser.add_argument('-c', '--create-datasets',
                        action='store_true',
                        help='Create new test datasets before testing')
    
    parser.add_argument('-r', '--repeats',
                        type=int,
                        default=3,
                        help='Number of times to repeat each test')
    
    parser.add_argument('-o', '--output',
                        default='benchmark_results.json',
                        help='Path to save benchmark results as JSON')
    
    parser.add_argument('-v', '--verbose',
                        action='store_true',
                        help='Enable verbose output')
    
    parser.add_argument('-s', '--sizes',
                        nargs='+',
                        type=int,
                        help='Specific dataset sizes to test')
    
    parser.add_argument('-p', '--paired-percents',
                        nargs='+',
                        type=float,
                        help='Specific paired percentages to test')
    
    parser.add_argument('-i', '--implementations',
                        nargs='+',
                        choices=['original', 'v2', 'v3'],
                        help='Specific implementations to test')
    
    return parser.parse_args()


def main() -> int:
    """Main function."""
    # Parse arguments
    args = parse_args()
    
    # Setup logger
    logger = get_logger("test_get_pairs")
    if args.verbose:
        logger.set_level(logging.DEBUG)
    
    try:
        # Install dependencies if needed
        try:
            import psutil
        except ImportError:
            logger.info("Installing psutil for memory tracking...")
            subprocess.check_call([sys.executable, "-m", "pip", "install", "psutil"])
            import psutil
            logger.info("psutil installed successfully")
        
        # Ensure the original script exists as needed
        if not os.path.exists(ORIGINAL_SCRIPT) and os.path.exists(V2_SCRIPT):
            logger.info(f"Backing up current implementation as {ORIGINAL_SCRIPT}")
            shutil.copy2(V2_SCRIPT, ORIGINAL_SCRIPT)
        
        # Use specified sizes and paired percentages, or defaults
        test_sizes = args.sizes if args.sizes else DEFAULT_TEST_SIZES
        paired_percents = args.paired_percents if args.paired_percents else DEFAULT_PAIRED_PERCENTS
        
        # Create test datasets if requested
        test_cases = {}
        if args.create_datasets:
            test_cases = create_test_datasets(args.test_dir, test_sizes, paired_percents)
        else:
            # Find existing test datasets
            test_dir_path = Path(args.test_dir)
            for size in test_sizes:
                for percent in paired_percents:
                    case_dir = test_dir_path / f"size_{size}_paired_{int(percent)}"
                    left_file = case_dir / f"reads_1_{size}_{int(percent)}.fastq"
                    right_file = case_dir / f"reads_2_{size}_{int(percent)}.fastq"
                    
                    if left_file.exists() and right_file.exists():
                        test_cases[f"{size}_{percent}"] = {
                            "left": str(left_file),
                            "right": str(right_file),
                            "size": size,
                            "paired_percent": percent
                        }
        
        if not test_cases:
            logger.error("No test datasets found. Use --create-datasets to create them.")
            return 1
        
        # Run benchmarks
        results = test_get_pairs(
            test_cases, 
            args.repeats, 
            args.output,
            args.implementations
        )
        
        # Analyze results
        analyze_results(results)
        
        return 0
    
    except Exception as e:
        logger.error(f"Error during testing: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(main()) 