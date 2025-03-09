"""Unit tests for the fastq_utils module."""

import os
import pytest
from src.utils.fastq_utils import (
    read_fastq_entry,
    write_fastq_entry,
    get_read_name,
    is_paired_entry,
    is_paired,
    open_file
)


def test_get_read_name():
    """Test the get_read_name function."""
    # Test normal read names
    assert get_read_name("@read1 description") == "read1"
    assert get_read_name("@read2/1 description") == "read2"
    assert get_read_name("@SRR123456.1 description") == "SRR123456.1"
    
    # Test with spaces in description
    assert get_read_name("@read3 this is a description") == "read3"
    
    # Test with no description
    assert get_read_name("@read4") == "read4"


def test_is_paired_entry():
    """Test the is_paired_entry function."""
    # Test paired entries
    assert is_paired_entry("@read1/1", "@read1/2") is True
    assert is_paired_entry("@SRR123.1 1", "@SRR123.1 2") is True
    
    # Test unpaired entries
    assert is_paired_entry("@read1/1", "@read2/2") is False
    assert is_paired_entry("@read1", "@read1/2") is False
    assert is_paired_entry("@read1/1", "@read1") is False


def test_is_paired():
    """Test the is_paired function."""
    # Test paired read names
    assert is_paired("read1/1", "read1/2") is True
    assert is_paired("SRR123.1", "SRR123.1") is True
    
    # Test unpaired read names
    assert is_paired("read1", "read2") is False
    assert is_paired("read1/1", "read2/2") is False


def test_read_fastq_entry(sample_fastq_file):
    """Test the read_fastq_entry function."""
    with open(sample_fastq_file, "r") as f:
        # Read the first entry
        entry = read_fastq_entry(f)
        assert entry[0] == "@seq1 description 1"
        assert entry[1] == "ACGTACGTACGT"
        assert entry[2] == "+"
        assert entry[3] == "IIIIIIIIIIII"
        
        # Read the second entry
        entry = read_fastq_entry(f)
        assert entry[0] == "@seq2 description 2"
        assert entry[1] == "TGCATGCATGCA"
        assert entry[2] == "+"
        assert entry[3] == "IIIIIIIIIIII"
        
        # Read the third entry
        entry = read_fastq_entry(f)
        assert entry[0] == "@seq3 description 3"
        assert entry[1] == "AAACCCTTTGGG"
        assert entry[2] == "+"
        assert entry[3] == "IIIIIIIIIIII"
        
        # There should be no more entries
        entry = read_fastq_entry(f)
        assert entry is None


def test_write_fastq_entry(temp_dir):
    """Test the write_fastq_entry function."""
    # Create a test file
    output_file = os.path.join(temp_dir, "output.fastq")
    with open(output_file, "w") as f:
        # Write a FASTQ entry
        write_fastq_entry(f, ["@test1", "ACGT", "+", "IIII"])
        write_fastq_entry(f, ["@test2", "TGCA", "+", "IIII"])
    
    # Verify the file content
    with open(output_file, "r") as f:
        content = f.read()
        assert content == "@test1\nACGT\n+\nIIII\n@test2\nTGCA\n+\nIIII\n"


def test_open_file(temp_dir):
    """Test the open_file function."""
    # Test opening a normal file
    file_path = os.path.join(temp_dir, "test.txt")
    with open(file_path, "w") as f:
        f.write("test content")
    
    with open_file(file_path) as f:
        content = f.read()
        assert content == "test content"
    
    # Test that the file is properly closed
    assert f.closed 