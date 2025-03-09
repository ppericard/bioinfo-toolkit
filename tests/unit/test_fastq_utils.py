"""Unit tests for the fastq_utils module."""

import os
import pytest
from src.utils.fastq_utils import (
    read_fastq_records,
    extract_read_id,
    open_file
)


def test_extract_read_id():
    """Test the extract_read_id function."""
    # Test normal read names
    assert extract_read_id("@read1 description") == "read1"
    assert extract_read_id("@read2/1 description") == "read2"
    assert extract_read_id("@SRR123456.1 description") == "SRR123456.1"
    
    # Test with spaces in description
    assert extract_read_id("@read3 this is a description") == "read3"
    
    # Test with no description
    assert extract_read_id("@read4") == "read4"


def test_read_fastq_records(sample_fastq_file):
    """Test the read_fastq_records function."""
    records = list(read_fastq_records(sample_fastq_file))
    
    # Check that we got the expected number of records
    assert len(records) == 3
    
    # Check the first record
    record_id, lines = records[0]
    assert record_id == "seq1"  # The ID part only
    assert lines[0] == "@seq1 description 1"  # Full header line
    assert lines[1] == "ACGTACGTACGT"
    assert lines[2] == "+"
    assert lines[3] == "IIIIIIIIIIII"
    
    # Check the second record
    record_id, lines = records[1]
    assert record_id == "seq2"  # The ID part only
    assert lines[0] == "@seq2 description 2"  # Full header line
    assert lines[1] == "TGCATGCATGCA"
    assert lines[2] == "+"
    assert lines[3] == "IIIIIIIIIIII"
    
    # Check the third record
    record_id, lines = records[2]
    assert record_id == "seq3"  # The ID part only
    assert lines[0] == "@seq3 description 3"  # Full header line
    assert lines[1] == "AAACCCTTTGGG"
    assert lines[2] == "+"
    assert lines[3] == "IIIIIIIIIIII"


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