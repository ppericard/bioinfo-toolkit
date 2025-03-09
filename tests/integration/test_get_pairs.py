"""Integration tests for the get_pairs script."""

import os
import pytest
import subprocess
import sys
from pathlib import Path

from bioinfotoolkit.scripts.fastq.get_pairs import get_pairs


def test_get_pairs_function(paired_fastq_files, temp_dir):
    """Test the get_pairs function directly."""
    left_file, right_file = paired_fastq_files
    output_dir = os.path.join(temp_dir, "output")
    os.makedirs(output_dir, exist_ok=True)
    
    # Run the get_pairs function
    stats = get_pairs(left_file, right_file, output_dir)
    
    # Check the output files
    assert os.path.exists(os.path.join(output_dir, "left.paired.fastq"))
    assert os.path.exists(os.path.join(output_dir, "right.paired.fastq"))
    assert os.path.exists(os.path.join(output_dir, "left.unpaired.fastq"))
    assert os.path.exists(os.path.join(output_dir, "right.unpaired.fastq"))
    
    # Check the stats
    assert stats["left"]["total"] == 3  # 3 reads in total
    assert stats["right"]["total"] == 3  # 3 reads in total
    assert stats["left"]["paired"] == 2  # 2 paired reads (read1 and read2)
    assert stats["right"]["paired"] == 2  # 2 paired reads (read1 and read2)
    assert stats["left"]["unpaired"] == 1  # 1 unpaired read (read3)
    assert stats["right"]["unpaired"] == 1  # 1 unpaired read (read4)


def test_get_pairs_script_execution(paired_fastq_files, temp_dir):
    """Test the get_pairs script execution through subprocess."""
    left_file, right_file = paired_fastq_files
    output_dir = os.path.join(temp_dir, "script_output")
    os.makedirs(output_dir, exist_ok=True)
    
    # Get the script path
    script_path = Path(__file__).parent.parent.parent / "src" / "scripts" / "fastq" / "get_pairs.py"
    
    # Run the script
    command = [
        sys.executable,
        str(script_path),
        "-l", left_file,
        "-r", right_file,
        "-o", output_dir
    ]
    
    result = subprocess.run(command, capture_output=True, text=True)
    
    # Check the process exit code
    assert result.returncode == 0
    
    # Check the output files
    assert os.path.exists(os.path.join(output_dir, "left.paired.fastq"))
    assert os.path.exists(os.path.join(output_dir, "right.paired.fastq"))
    assert os.path.exists(os.path.join(output_dir, "left.unpaired.fastq"))
    assert os.path.exists(os.path.join(output_dir, "right.unpaired.fastq"))
    
    # Check the output for expected messages
    # The output format is "Left file: X paired, Y unpaired"
    # The logging output goes to stderr, not stdout
    assert "Left file: 2 paired, 1 unpaired" in result.stderr
    assert "Right file: 2 paired, 1 unpaired" in result.stderr 