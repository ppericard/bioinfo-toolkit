# Migration Plan for Bioinfo-Toolkit

This document outlines the steps needed to complete the reorganization of the bioinfo-toolkit project.

## Changes Already Made

1. Created a proper Python package structure with `pyproject.toml`
2. Set up the main package directory at `src/bioinfotoolkit/`
3. Created a command-line interface module at `src/bioinfotoolkit/cli.py`
4. Created a simplified entry point script at `bioinfo-toolkit`
5. Updated `.gitignore` to include additional patterns
6. Created a simplified README that reflects the new structure
7. Tested the basic functionality of the new structure

## Remaining Tasks

### 1. Move Script Files

Move all script files from the old structure to the new structure:

```bash
# For each script category (fasta, fastq, conversion, tools, benchmark)
cp src/scripts/fasta/*.py src/bioinfotoolkit/scripts/fasta/
cp src/scripts/fastq/*.py src/bioinfotoolkit/scripts/fastq/
cp src/scripts/conversion/*.py src/bioinfotoolkit/scripts/conversion/
cp src/scripts/tools/*.py src/bioinfotoolkit/scripts/tools/
cp src/scripts/benchmark/*.py src/bioinfotoolkit/scripts/benchmark/
```

### 2. Move Utility Modules

Move all utility modules from the old structure to the new structure:

```bash
cp src/utils/*.py src/bioinfotoolkit/utils/
```

### 3. Update Imports in Script Files

Update import statements in all script files to use the new package structure. For example:

- Change `from utils.bioinfo_logger import get_logger` to `from bioinfotoolkit.utils.bioinfo_logger import get_logger`
- Change `from utils.fastq_utils import read_fastq` to `from bioinfotoolkit.utils.fastq_utils import read_fastq`

### 4. Update the Utils __init__.py File

Update the `src/bioinfotoolkit/utils/__init__.py` file to import all utility functions and classes:

```python
"""Utility modules for bioinfo-toolkit."""

from bioinfotoolkit.utils.bioinfo_logger import get_logger, BioinfLogger
from bioinfotoolkit.utils.fastq_utils import read_fastq, write_fastq
from bioinfotoolkit.utils.memory_tracker import MemoryTracker

__all__ = [
    'get_logger', 'BioinfLogger',
    'read_fastq', 'write_fastq',
    'MemoryTracker'
]
```

### 5. Test All Scripts

Test all scripts to ensure they work with the new package structure:

```bash
bioinfo-toolkit fasta_length_filter --help
bioinfo-toolkit fastq_to_fasta --help
# Test other scripts as needed
```

### 6. Clean Up Old Files

Once everything is working with the new structure, remove the old files:

```bash
# Remove old script files
rm -rf src/scripts/

# Remove old utility modules
rm -rf src/utils/

# Remove old entry point script
rm bioinfo-toolkit.py

# Remove old requirements files
rm requirements.txt requirements-dev.txt
```

### 7. Update Tests

Update import statements in test files to use the new package structure.

### 8. Commit Changes

Commit the changes in logical units:

```bash
# Commit the new package structure
git add pyproject.toml src/bioinfotoolkit/ bioinfo-toolkit
git commit -m "Implement proper Python package structure"

# Commit the updated scripts
git add src/bioinfotoolkit/scripts/
git commit -m "Move scripts to new package structure"

# Commit the updated utility modules
git add src/bioinfotoolkit/utils/
git commit -m "Move utility modules to new package structure"

# Commit the cleanup
git add -u
git commit -m "Remove old files after migration to new package structure"
```

## Benefits of the New Structure

1. **Standard Python Package**: The project now follows standard Python packaging conventions
2. **Easier Installation**: Users can install the package with pip
3. **Better Imports**: Cleaner import statements with a consistent package namespace
4. **Simplified Entry Point**: A single command-line entry point for all scripts
5. **Better Dependency Management**: Dependencies are managed in pyproject.toml
6. **Development Mode**: Easy development with pip's editable mode 