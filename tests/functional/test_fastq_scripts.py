"""Functional tests for FASTQ processing scripts."""

import os
import pytest
import subprocess
import sys
from pathlib import Path


def test_fastq_name_filter(sample_fastq_file, temp_dir):
    """Test the fastq_name_filter script."""
    # Get script path
    script_path = Path(__file__).parent.parent.parent / "src" / "scripts" / "fastq" / "fastq_name_filter.py"
    
    # Create a patterns file
    patterns_file = os.path.join(temp_dir, "patterns.txt")
    with open(patterns_file, "w") as f:
        f.write("seq1\nseq3\n")
    
    # Run the script
    output_file = os.path.join(temp_dir, "filtered_names.fastq")
    cmd = [
        sys.executable,
        str(script_path),
        "-i", sample_fastq_file,
        "-o", output_file,
        "-f", patterns_file
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    # Check that the script executed successfully
    assert result.returncode == 0
    assert os.path.exists(output_file)
    
    # Check the output file content
    with open(output_file, "r") as f:
        content = f.read()
        # Only seq1 and seq3 should be in the output
        assert "@seq1" in content
        assert "@seq3" in content
        assert "@seq2" not in content


def test_fastq_umi_merge(temp_dir):
    """Test the fastq_umi_merge script."""
    # Create test files with UMIs
    umi_file = os.path.join(temp_dir, "umis.fastq")
    with open(umi_file, "w") as f:
        f.write("@read1\nAAAA\n+\nIIII\n")
        f.write("@read2\nCCCC\n+\nIIII\n")
    
    reads_file = os.path.join(temp_dir, "reads.fastq")
    with open(reads_file, "w") as f:
        f.write("@read1\nACGTACGTACGT\n+\nIIIIIIIIIIII\n")
        f.write("@read2\nTGCATGCATGCA\n+\nIIIIIIIIIIII\n")
    
    # Get script path
    script_path = Path(__file__).parent.parent.parent / "src" / "scripts" / "fastq" / "fastq_umi_merge.py"
    
    # Run the script
    output_file = os.path.join(temp_dir, "merged.fastq")
    cmd = [
        sys.executable,
        str(script_path),
        "-u", umi_file,
        "-i", reads_file,
        "-o", output_file
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    # Check that the script executed successfully
    assert result.returncode == 0
    assert os.path.exists(output_file)
    
    # Check the output file content
    with open(output_file, "r") as f:
        content = f.read()
        # UMI should be prepended to the read sequence
        assert "AAAAACGTACGTACGT" in content
        assert "CCCCTGCATGCATGCA" in content


def test_gener_sample_fastq(sample_fastq_file, temp_dir):
    """Test the gener_sample_fastq script."""
    # Get script path
    script_path = Path(__file__).parent.parent.parent / "src" / "scripts" / "fastq" / "gener_sample_fastq.py"
    
    # Run the script
    output_file = os.path.join(temp_dir, "sample.fastq")
    cmd = [
        sys.executable,
        str(script_path),
        "-i", sample_fastq_file,
        "-o", output_file,
        "-p", "1.0"  # Use 100% probability to ensure we get all sequences
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    # Check that the script executed successfully
    assert result.returncode == 0
    assert os.path.exists(output_file)
    
    # Check the output file content
    with open(output_file, "r") as f:
        content = f.read()
        # Should contain at least one sequence
        assert "@" in content


def test_split_paired_fastq(temp_dir):
    """Test the split_paired_fastq script."""
    # Create a paired FASTQ file (interleaved format)
    input_file = os.path.join(temp_dir, "interleaved.fastq")
    with open(input_file, "w") as f:
        # First pair
        f.write("@read1/1\nACGTACGTACGT\n+\nIIIIIIIIIIII\n")
        f.write("@read1/2\nTGCATGCATGCA\n+\nIIIIIIIIIIII\n")
        # Second pair
        f.write("@read2/1\nGTACGTACGTAC\n+\nIIIIIIIIIIII\n")
        f.write("@read2/2\nCATGCATGCATG\n+\nIIIIIIIIIIII\n")
    
    # Get script path
    script_path = Path(__file__).parent.parent.parent / "src" / "scripts" / "fastq" / "split_paired_fastq.py"
    
    # Run the script
    output_left = os.path.join(temp_dir, "split_1.fastq")
    output_right = os.path.join(temp_dir, "split_2.fastq")
    cmd = [
        sys.executable,
        str(script_path),
        "-i", input_file,
        "-1", output_left,
        "-2", output_right
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    # Check that the script executed successfully
    assert result.returncode == 0
    
    # Check that the output files exist and contain the correct reads
    assert os.path.exists(output_left)
    assert os.path.exists(output_right)
    
    with open(output_left, "r") as f:
        content = f.read()
        assert "@read1/1" in content
        assert "@read2/1" in content
        assert "@read1/2" not in content
        assert "@read2/2" not in content
    
    with open(output_right, "r") as f:
        content = f.read()
        assert "@read1/2" in content
        assert "@read2/2" in content
        assert "@read1/1" not in content
        assert "@read2/1" not in content 