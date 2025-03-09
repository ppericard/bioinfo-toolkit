"""Common fixtures and utilities for bioinfo-toolkit tests."""

import os
import pytest
import tempfile
import shutil
from pathlib import Path

@pytest.fixture
def temp_dir():
    """Create a temporary directory for test files."""
    temp_dir = tempfile.mkdtemp()
    yield temp_dir
    shutil.rmtree(temp_dir)

@pytest.fixture
def data_dir():
    """Return the path to the test data directory."""
    base_dir = Path(__file__).parent
    return os.path.join(base_dir, "data")

@pytest.fixture
def sample_fasta_file(temp_dir):
    """Create a sample FASTA file for testing."""
    fasta_content = (
        ">seq1 description 1\n"
        "ACGTACGTACGT\n"
        ">seq2 description 2\n"
        "TGCATGCATGCA\n"
        ">seq3 description 3\n"
        "AAACCCTTTGGG\n"
    )
    file_path = os.path.join(temp_dir, "sample.fasta")
    with open(file_path, "w") as f:
        f.write(fasta_content)
    return file_path

@pytest.fixture
def sample_fastq_file(temp_dir):
    """Create a sample FASTQ file for testing."""
    fastq_content = (
        "@seq1 description 1\n"
        "ACGTACGTACGT\n"
        "+\n"
        "IIIIIIIIIIII\n"
        "@seq2 description 2\n"
        "TGCATGCATGCA\n"
        "+\n"
        "IIIIIIIIIIII\n"
        "@seq3 description 3\n"
        "AAACCCTTTGGG\n"
        "+\n"
        "IIIIIIIIIIII\n"
    )
    file_path = os.path.join(temp_dir, "sample.fastq")
    with open(file_path, "w") as f:
        f.write(fastq_content)
    return file_path

@pytest.fixture
def paired_fastq_files(temp_dir):
    """Create paired FASTQ files for testing."""
    # Create left reads
    left_content = (
        "@read1/1\n"
        "ACGTACGTACGT\n"
        "+\n"
        "IIIIIIIIIIII\n"
        "@read2/1\n"
        "TGCATGCATGCA\n"
        "+\n"
        "IIIIIIIIIIII\n"
        "@read3/1\n"
        "AAACCCTTTGGG\n"
        "+\n"
        "IIIIIIIIIIII\n"
    )
    left_path = os.path.join(temp_dir, "sample_1.fastq")
    with open(left_path, "w") as f:
        f.write(left_content)
    
    # Create right reads
    right_content = (
        "@read1/2\n"
        "TGCATGCATGCA\n"
        "+\n"
        "IIIIIIIIIIII\n"
        "@read2/2\n"
        "ACGTACGTACGT\n"
        "+\n"
        "IIIIIIIIIIII\n"
        "@read4/2\n"  # Unpaired read
        "CCCGGGAAATTT\n"
        "+\n"
        "IIIIIIIIIIII\n"
    )
    right_path = os.path.join(temp_dir, "sample_2.fastq")
    with open(right_path, "w") as f:
        f.write(right_content)
    
    return left_path, right_path 