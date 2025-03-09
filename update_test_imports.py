#!/usr/bin/env python3
"""
Script to update import statements in all test files.
This will change imports from:
- from src.utils.x import y
to:
- from bioinfotoolkit.utils.x import y

And from:
- import src.utils.x
to:
- import bioinfotoolkit.utils.x

And from:
- from src.scripts.x import y
to:
- from bioinfotoolkit.scripts.x import y
"""

import os
import re
from pathlib import Path

def update_file_imports(file_path):
    """Update imports in a single file."""
    with open(file_path, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # Update imports
    # Pattern 1: from src.utils.x import y
    pattern1 = r'from\s+src\.utils\.(.*?)\s+import\s+(.*?)$'
    replacement1 = r'from bioinfotoolkit.utils.\1 import \2'
    content = re.sub(pattern1, replacement1, content, flags=re.MULTILINE)
    
    # Pattern 2: import src.utils.x
    pattern2 = r'import\s+src\.utils\.(.*?)$'
    replacement2 = r'import bioinfotoolkit.utils.\1'
    content = re.sub(pattern2, replacement2, content, flags=re.MULTILINE)

    # Pattern 3: from src.scripts.x import y
    pattern3 = r'from\s+src\.scripts\.(.*?)\s+import\s+(.*?)$'
    replacement3 = r'from bioinfotoolkit.scripts.\1 import \2'
    content = re.sub(pattern3, replacement3, content, flags=re.MULTILINE)
    
    # Pattern 4: import src.scripts.x
    pattern4 = r'import\s+src\.scripts\.(.*?)$'
    replacement4 = r'import bioinfotoolkit.scripts.\1'
    content = re.sub(pattern4, replacement4, content, flags=re.MULTILINE)
    
    # Write updated content back to file
    with open(file_path, 'w', encoding='utf-8') as f:
        f.write(content)
    
    print(f"Updated imports in {file_path}")

def main():
    """Main function to update imports in all test files."""
    base_dir = Path('tests')
    
    # Process each test directory
    for test_dir in ['unit', 'integration', 'functional', 'performance']:
        test_path = base_dir / test_dir
        if not test_path.exists():
            continue
        
        for test_file in test_path.glob('*.py'):
            if test_file.name == '__init__.py':
                continue
            update_file_imports(test_file)
    
    # Also update conftest.py
    conftest_path = base_dir / 'conftest.py'
    if conftest_path.exists():
        update_file_imports(conftest_path)
    
    print("All test imports updated successfully!")

if __name__ == "__main__":
    main() 