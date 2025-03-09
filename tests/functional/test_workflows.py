"""Functional tests for common bioinformatics workflows."""

import os
import pytest
import subprocess
import sys
from pathlib import Path


def test_fastq_to_fasta_filter_workflow(sample_fastq_file, temp_dir):
    """Test a workflow: convert FASTQ to FASTA and then filter by length."""
    # Get script paths
    fastq_to_fasta_path = Path(__file__).parent.parent.parent / "src" / "scripts" / "conversion" / "fastq_to_fasta.py"
    fasta_filter_path = Path(__file__).parent.parent.parent / "src" / "scripts" / "fasta" / "fasta_length_filter.py"
    
    # Define output files
    fasta_output = os.path.join(temp_dir, "output.fasta")
    filtered_output = os.path.join(temp_dir, "filtered.fasta")
    
    # Step 1: Convert FASTQ to FASTA
    fastq_to_fasta_cmd = [
        sys.executable,
        str(fastq_to_fasta_path),
        "-i", sample_fastq_file,
        "-o", fasta_output
    ]
    
    result1 = subprocess.run(fastq_to_fasta_cmd, capture_output=True, text=True)
    assert result1.returncode == 0
    assert os.path.exists(fasta_output)
    
    # Step 2: Filter FASTA by length
    fasta_filter_cmd = [
        sys.executable,
        str(fasta_filter_path),
        "-i", fasta_output,
        "-o", filtered_output,
        "-m", "10",  # Minimum length
        "-M", "20"   # Maximum length
    ]
    
    result2 = subprocess.run(fasta_filter_cmd, capture_output=True, text=True)
    assert result2.returncode == 0
    assert os.path.exists(filtered_output)
    
    # Verify the content of the filtered file
    with open(filtered_output, "r") as f:
        content = f.read()
        # All sequences are 12 nucleotides, so all should pass the filter
        assert ">seq1" in content
        assert ">seq2" in content
        assert ">seq3" in content


def test_paired_reads_workflow(paired_fastq_files, temp_dir):
    """Test a workflow: get paired reads and convert to FASTA."""
    # Get script paths
    get_pairs_path = Path(__file__).parent.parent.parent / "src" / "scripts" / "fastq" / "get_pairs.py"
    fastq_to_fasta_path = Path(__file__).parent.parent.parent / "src" / "scripts" / "conversion" / "fastq_to_fasta.py"
    
    # Define output directories and files
    pairs_output_dir = os.path.join(temp_dir, "pairs_output")
    os.makedirs(pairs_output_dir, exist_ok=True)
    fasta_output = os.path.join(temp_dir, "paired.fasta")
    
    # Step 1: Get paired reads
    left_file, right_file = paired_fastq_files
    get_pairs_cmd = [
        sys.executable,
        str(get_pairs_path),
        "-l", left_file,
        "-r", right_file,
        "-o", pairs_output_dir
    ]
    
    result1 = subprocess.run(get_pairs_cmd, capture_output=True, text=True)
    assert result1.returncode == 0
    
    left_paired_file = os.path.join(pairs_output_dir, "left.paired.fastq")
    assert os.path.exists(left_paired_file)
    
    # Step 2: Convert paired reads to FASTA
    fastq_to_fasta_cmd = [
        sys.executable,
        str(fastq_to_fasta_path),
        "-i", left_paired_file,
        "-o", fasta_output
    ]
    
    result2 = subprocess.run(fastq_to_fasta_cmd, capture_output=True, text=True)
    assert result2.returncode == 0
    assert os.path.exists(fasta_output)
    
    # Verify the content of the FASTA file
    with open(fasta_output, "r") as f:
        content = f.read()
        # Check that the paired reads are present
        assert ">read1" in content
        assert ">read2" in content
        # read3 should not be present as it was unpaired
        assert ">read3" not in content 