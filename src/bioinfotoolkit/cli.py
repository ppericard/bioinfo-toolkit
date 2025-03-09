#!/usr/bin/env python3
"""
Command-line interface for bioinfo-toolkit.

This module serves as an entry point to run the various tools in the toolkit.
"""

import sys
import argparse
import importlib.util
from pathlib import Path
from typing import Dict, List, Optional

from bioinfotoolkit import __version__

# Define script categories and their descriptions
SCRIPT_CATEGORIES = {
    'fasta': 'FASTA processing tools',
    'fastq': 'FASTQ processing tools',
    'conversion': 'File format conversion tools',
    'tools': 'Utility tools'
}

def get_script_paths() -> Dict[str, Path]:
    """
    Get a mapping of script names to their paths.
    
    Returns:
        Dictionary mapping script names to their file paths
    """
    scripts = {}
    package_path = Path(__file__).parent / "scripts"
    
    for category in SCRIPT_CATEGORIES:
        category_dir = package_path / category
        if not category_dir.exists():
            continue
        
        for script_file in category_dir.glob("*.py"):
            if script_file.stem.startswith('__'):
                continue
            
            scripts[script_file.stem] = script_file
    
    return scripts

def list_available_scripts() -> None:
    """List all available scripts in the toolkit, organized by category."""
    package_path = Path(__file__).parent / "scripts"
    
    print("Available scripts in bioinfo-toolkit:")
    print("=====================================")
    
    for category, description in SCRIPT_CATEGORIES.items():
        print(f"\n{description}:")
        print("-" * len(description))
        
        category_dir = package_path / category
        if not category_dir.exists():
            print(f"  No scripts found in {category}")
            continue
            
        scripts = [f.stem for f in category_dir.glob("*.py") 
                  if not f.stem.startswith('__')]
        
        if not scripts:
            print("  No scripts found")
            continue
        
        for script in sorted(scripts):
            print(f"  {script}")

def execute_script(script_name: str, script_args: List[str]) -> int:
    """
    Execute a script with the given arguments.
    
    Args:
        script_name: Name of the script to execute
        script_args: Arguments to pass to the script
        
    Returns:
        Exit code of the script
    """
    scripts = get_script_paths()
    
    if script_name not in scripts:
        print(f"Error: Script '{script_name}' not found.")
        list_available_scripts()
        return 1
    
    script_path = scripts[script_name]
    
    # Set up the arguments for the script
    sys.argv = [str(script_path)] + script_args
    
    # Load and execute the script
    spec = importlib.util.spec_from_file_location(script_name, script_path)
    if spec is None:
        print(f"Error: Could not load script '{script_name}'.")
        return 1
    
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    
    # Call the main function of the script
    if hasattr(module, 'main'):
        return module.main()
    
    return 0

def main() -> int:
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
    
    return execute_script(args.script, args.args)

if __name__ == "__main__":
    sys.exit(main()) 