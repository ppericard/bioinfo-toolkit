#!/usr/bin/env python3
"""
Script to update import statements in all script files.
This will change imports from:
- from utils.x import y
to:
- from bioinfotoolkit.utils.x import y

And from:
- import utils.x
to:
- import bioinfotoolkit.utils.x
"""

import os
import re
from pathlib import Path

def update_file_imports(file_path):
    """Update imports in a single file."""
    with open(file_path, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # Update imports
    # Pattern 1: from utils.x import y
    pattern1 = r'from\s+utils\.(.*?)\s+import\s+(.*?)$'
    replacement1 = r'from bioinfotoolkit.utils.\1 import \2'
    content = re.sub(pattern1, replacement1, content, flags=re.MULTILINE)
    
    # Pattern 2: import utils.x
    pattern2 = r'import\s+utils\.(.*?)$'
    replacement2 = r'import bioinfotoolkit.utils.\1'
    content = re.sub(pattern2, replacement2, content, flags=re.MULTILINE)

    # Pattern 3: from src.utils.x import y
    pattern3 = r'from\s+src\.utils\.(.*?)\s+import\s+(.*?)$'
    replacement3 = r'from bioinfotoolkit.utils.\1 import \2'
    content = re.sub(pattern3, replacement3, content, flags=re.MULTILINE)
    
    # Write updated content back to file
    with open(file_path, 'w', encoding='utf-8') as f:
        f.write(content)
    
    print(f"Updated imports in {file_path}")

def main():
    """Main function to update imports in all script files."""
    base_dir = Path('src/bioinfotoolkit/scripts')
    
    for category in ['fasta', 'fastq', 'conversion', 'tools', 'benchmark']:
        category_dir = base_dir / category
        if not category_dir.exists():
            continue
        
        for script_file in category_dir.glob('*.py'):
            if script_file.name == '__init__.py':
                continue
            update_file_imports(script_file)
    
    # Also update utility modules
    utils_dir = Path('src/bioinfotoolkit/utils')
    for util_file in utils_dir.glob('*.py'):
        if util_file.name == '__init__.py':
            continue
        update_file_imports(util_file)
    
    print("All imports updated successfully!")

if __name__ == "__main__":
    main() 