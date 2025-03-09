#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Benchmark Runner for bioinfo-toolkit

This script runs all benchmarks and generates reports.
"""

import os
import sys
import argparse
import importlib
import subprocess
from pathlib import Path

def list_benchmarks():
    """List all available benchmarks in the benchmarks directory."""
    benchmark_dir = Path(__file__).parent
    benchmarks = []
    
    for file in benchmark_dir.glob("*.py"):
        if file.name == "__init__.py" or file.name == "run_benchmarks.py":
            continue
        
        name = file.stem
        if name.startswith("benchmark_"):
            benchmarks.append(name)
    
    return benchmarks

def run_benchmark(benchmark_name, args=None):
    """Run a specific benchmark."""
    if args is None:
        args = []
    
    benchmark_path = Path(__file__).parent / f"{benchmark_name}.py"
    
    if not benchmark_path.exists():
        print(f"Benchmark '{benchmark_name}' not found.")
        return 1
    
    cmd = [sys.executable, str(benchmark_path)] + args
    return subprocess.call(cmd)

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Run benchmarks for bioinfo-toolkit",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        'benchmark', 
        nargs='?',
        help='Benchmark to run. Omit to list available benchmarks.'
    )
    
    parser.add_argument(
        'args', 
        nargs=argparse.REMAINDER,
        help='Arguments to pass to the benchmark.'
    )
    
    parser.add_argument(
        '--all', 
        action='store_true',
        help='Run all benchmarks.'
    )
    
    return parser.parse_args()

def main():
    """Main function."""
    args = parse_args()
    
    if args.all:
        benchmarks = list_benchmarks()
        print(f"Running all {len(benchmarks)} benchmarks...")
        
        for benchmark in benchmarks:
            print(f"\nRunning benchmark: {benchmark}")
            run_benchmark(benchmark)
        
        return 0
    
    if not args.benchmark:
        print("Available benchmarks:")
        for benchmark in list_benchmarks():
            print(f"  {benchmark}")
        return 0
    
    return run_benchmark(args.benchmark, args.args)

if __name__ == "__main__":
    sys.exit(main()) 