#!/usr/bin/env python3
"""
Command-line interface for bioinfo-toolkit.

This module serves as an entry point to run the various tools in the toolkit.
"""

import os
import sys
import argparse
import importlib.util
from pathlib import Path

from bioinfotoolkit import __version__

def list_available_scripts():
    """List all available scripts in the toolkit."""
    script_categories = {
        'fasta': 'FASTA processing tools',
        'fastq': 'FASTQ processing tools',
        'conversion': 'File format conversion tools',
        'tools': 'Utility tools'
    }
    
    print("Available scripts in bioinfo-toolkit:")
    print("=====================================")
    
    # Get the base directory where the modules are installed
    package_path = Path(__file__).parent / "scripts"
    
    for category, description in script_categories.items():
        print(f"\n{description}:")
        print("-" * len(description))
        
        category_dir = package_path / category
        if not category_dir.exists():
            print(f"  No scripts found in {category}")
            continue
            
        scripts = [f.stem for f in category_dir.glob("*.py") 
                  if not f.stem.startswith('__')]
        
        for script in sorted(scripts):
            print(f"  {script}")

def main():
    """Main entry point for the toolkit."""
    parser = argparse.ArgumentParser(
        description="Bioinfo-toolkit - A collection of Python scripts for bioinformatics analysis"
    )
    parser.add_argument(
        '--version', action='version', version=f'%(prog)s {__version__}'
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
    script_categories = ['fasta', 'fastq', 'conversion', 'tools']
    script_path = None
    
    package_path = Path(__file__).parent / "scripts"
    
    for category in script_categories:
        candidate_path = package_path / category / f"{args.script}.py"
        if candidate_path.exists():
            script_path = str(candidate_path)
            break
    
    if not script_path:
        print(f"Error: Script '{args.script}' not found in any category.")
        list_available_scripts()
        return 1
    
    # Run the script with the provided arguments
    sys.argv = [script_path] + args.args
    
    # Load and execute the script module
    spec = importlib.util.spec_from_file_location(args.script, script_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    
    # If the script has a main function, call it
    if hasattr(module, 'main'):
        return module.main()
    
    return 0

if __name__ == "__main__":
    sys.exit(main()) 