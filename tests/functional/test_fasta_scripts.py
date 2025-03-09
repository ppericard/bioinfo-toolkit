"""Functional tests for FASTA processing scripts."""

import os
import pytest
import subprocess
import sys
from pathlib import Path


def test_fasta_length_histo(sample_fasta_file, temp_dir):
    """Test the fasta_length_histo script."""
    # Get script path
    script_path = Path(__file__).parent.parent.parent / "src" / "scripts" / "fasta" / "fasta_length_histo.py"
    
    # Run the script
    output_file = os.path.join(temp_dir, "length_histo.pdf")
    cmd = [
        sys.executable,
        str(script_path),
        sample_fasta_file,  # Positional FASTA argument
        "-o", output_file
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    # Check that the script executed successfully
    assert result.returncode == 0
    assert os.path.exists(output_file)
    
    # Since we're generating a PDF, we can't easily check its content,
    # but at least we verified it was created


def test_fasta_n_filter(temp_dir):
    """Test the fasta_n_filter script."""
    # Create a test FASTA file with some N's
    input_file = os.path.join(temp_dir, "with_n.fasta")
    with open(input_file, "w") as f:
        f.write(">seq1\nACGTACGTACGT\n")
        f.write(">seq2\nACGTNNNNACGT\n")
        f.write(">seq3\nNNNTACGTACGT\n")
    
    # Get script path
    script_path = Path(__file__).parent.parent.parent / "src" / "scripts" / "fasta" / "fasta_n_filter.py"
    
    # Run the script
    output_file = os.path.join(temp_dir, "filtered_n.fasta")
    cmd = [
        sys.executable,
        str(script_path),
        "-i", input_file,
        "-o", output_file,
        "-r", "0.3"  # Max N rate of 30%
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    # Check that the script executed successfully
    assert result.returncode == 0
    assert os.path.exists(output_file)
    
    # Check the output file content
    with open(output_file, "r") as f:
        content = f.read()
        # seq1 has 0% N, seq2 has 33% N, seq3 has 25% N
        # With threshold at 30%, only seq1 and seq3 should pass
        assert ">seq1" in content
        assert ">seq3" in content
        assert ">seq2" not in content


def test_fasta_name_filter(sample_fasta_file, temp_dir):
    """Test the fasta_name_filter script."""
    # Get script path
    script_path = Path(__file__).parent.parent.parent / "src" / "scripts" / "fasta" / "fasta_name_filter.py"
    
    # Create a patterns file
    patterns_file = os.path.join(temp_dir, "patterns.txt")
    with open(patterns_file, "w") as f:
        f.write("seq1\nseq3\n")
    
    # Run the script
    output_file = os.path.join(temp_dir, "filtered_names.fasta")
    cmd = [
        sys.executable,
        str(script_path),
        "-i", sample_fasta_file,
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
        assert ">seq1" in content
        assert ">seq3" in content
        assert ">seq2" not in content


def test_gener_sample_fasta(sample_fasta_file, temp_dir):
    """Test the gener_sample_fasta script."""
    # Get script path
    script_path = Path(__file__).parent.parent.parent / "src" / "scripts" / "fasta" / "gener_sample_fasta.py"
    
    # Run the script
    output_file = os.path.join(temp_dir, "sample.fasta")
    cmd = [
        sys.executable,
        str(script_path),
        "-i", sample_fasta_file,
        "-o", output_file,
        "-p", "0.5"  # Use probability instead of fixed number
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    # Check that the script executed successfully
    assert result.returncode == 0
    assert os.path.exists(output_file)
    
    # Check the output file content
    with open(output_file, "r") as f:
        content = f.read()
        # Should contain at least one sequence
        assert ">" in content


def test_sort_fasta_by_length(sample_fasta_file, temp_dir):
    """Test the sort_fasta_by_length script."""
    # Get script path
    script_path = Path(__file__).parent.parent.parent / "src" / "scripts" / "fasta" / "sort_fasta_by_length.py"
    
    # Run the script
    output_file = os.path.join(temp_dir, "sorted.fasta")
    cmd = [
        sys.executable,
        str(script_path),
        "--input", sample_fasta_file,
        "--output", output_file
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    # Check that the script executed successfully
    assert result.returncode == 0
    assert os.path.exists(output_file)
    
    # Since all test sequences are the same length, just verify all sequences are present
    with open(output_file, "r") as f:
        content = f.read()
        assert ">seq1" in content
        assert ">seq2" in content
        assert ">seq3" in content 