#!/usr/bin/env python3
"""
Bioinfo-toolkit - A collection of Python scripts for bioinformatics analysis.

This script serves as an entry point to run the various tools in the toolkit.
"""

import os
import sys
import argparse
import importlib
import pkgutil

def list_available_scripts():
    """List all available scripts in the toolkit."""
    script_categories = {
        'fasta': 'FASTA processing tools',
        'fastq': 'FASTQ processing tools',
        'conversion': 'File format conversion tools',
        'tools': 'Utility tools',
        'benchmark': 'Benchmarking and testing tools'
    }
    
    print("Available scripts in bioinfo-toolkit:")
    print("=====================================")
    
    for category, description in script_categories.items():
        print(f"\n{description}:")
        print("-" * len(description))
        
        category_dir = os.path.join('src', 'scripts', category)
        if not os.path.exists(category_dir):
            print(f"  No scripts found in {category}")
            continue
            
        scripts = [f for f in os.listdir(category_dir) 
                  if f.endswith('.py') and not f.startswith('__')]
        
        for script in sorted(scripts):
            script_name = script[:-3]  # Remove .py extension
            print(f"  {script_name}")

def main():
    """Main entry point for the toolkit."""
    parser = argparse.ArgumentParser(
        description="Bioinfo-toolkit - A collection of Python scripts for bioinformatics analysis"
    )
    parser.add_argument(
        'script', nargs='?', help='Script to run (omit to list available scripts)'
    )
    parser.add_argument(
        'args', nargs=argparse.REMAINDER, help='Arguments to pass to the script'
    )
    
    args = parser.parse_args()
    
    if not args.script:
        list_available_scripts()
        return 0
    
    # Try to find the script in the different categories
    script_categories = ['fasta', 'fastq', 'conversion', 'tools', 'benchmark']
    script_path = None
    
    for category in script_categories:
        candidate_path = os.path.join('src', 'scripts', category, f"{args.script}.py")
        if os.path.exists(candidate_path):
            script_path = candidate_path
            break
    
    if not script_path:
        print(f"Error: Script '{args.script}' not found in any category.")
        list_available_scripts()
        return 1
    
    # Run the script with the provided arguments
    sys.argv = [script_path] + args.args
    with open(script_path, 'r') as f:
        exec(f.read())
    
    return 0

if __name__ == "__main__":
    sys.exit(main()) 