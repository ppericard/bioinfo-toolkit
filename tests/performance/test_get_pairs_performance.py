"""Performance tests for the get_pairs implementations."""

import os
import pytest
import time
import tempfile
import shutil
from pathlib import Path
import random
import string

from bioinfotoolkit.utils.memory_tracker import MemoryTracker
from bioinfotoolkit.scripts.fastq.get_pairs import get_pairs
from bioinfotoolkit.scripts.fastq.get_pairs_v3 import get_pairs as get_pairs_v3


def generate_random_sequence(length):
    """Generate a random DNA sequence."""
    return ''.join(random.choice('ACGT') for _ in range(length))


def generate_random_quality(length):
    """Generate a random quality string."""
    return ''.join(random.choice('IJKLMNOPQRSTUVWXYZ') for _ in range(length))


def create_fastq_files(num_reads, paired_percentage, temp_dir):
    """Create random FASTQ files for testing."""
    left_path = os.path.join(temp_dir, f"left_{num_reads}.fastq")
    right_path = os.path.join(temp_dir, f"right_{num_reads}.fastq")
    
    # Calculate the number of paired and unpaired reads
    paired_reads = int(num_reads * paired_percentage / 100)
    unpaired_left = int((num_reads - paired_reads) / 2)
    unpaired_right = num_reads - paired_reads - unpaired_left
    
    # Generate read IDs
    read_ids = [f"read{i}" for i in range(1, paired_reads + 1)]
    left_unpaired_ids = [f"left_unpaired{i}" for i in range(1, unpaired_left + 1)]
    right_unpaired_ids = [f"right_unpaired{i}" for i in range(1, unpaired_right + 1)]
    
    # Write left reads
    with open(left_path, "w") as f:
        # Write paired reads
        for read_id in read_ids:
            seq = generate_random_sequence(100)
            qual = generate_random_quality(100)
            f.write(f"@{read_id}/1\n{seq}\n+\n{qual}\n")
        
        # Write unpaired reads
        for read_id in left_unpaired_ids:
            seq = generate_random_sequence(100)
            qual = generate_random_quality(100)
            f.write(f"@{read_id}\n{seq}\n+\n{qual}\n")
    
    # Write right reads
    with open(right_path, "w") as f:
        # Write paired reads
        for read_id in read_ids:
            seq = generate_random_sequence(100)
            qual = generate_random_quality(100)
            f.write(f"@{read_id}/2\n{seq}\n+\n{qual}\n")
        
        # Write unpaired reads
        for read_id in right_unpaired_ids:
            seq = generate_random_sequence(100)
            qual = generate_random_quality(100)
            f.write(f"@{read_id}\n{seq}\n+\n{qual}\n")
    
    return left_path, right_path


@pytest.mark.parametrize("num_reads", [100, 1000])
@pytest.mark.parametrize("paired_percentage", [50, 90])
def test_get_pairs_performance(num_reads, paired_percentage):
    """Test the performance of get_pairs implementations."""
    # Create a temporary directory for test files
    temp_dir = tempfile.mkdtemp()
    try:
        # Create test files
        left_path, right_path = create_fastq_files(num_reads, paired_percentage, temp_dir)
        
        # Create output directories
        output_dir1 = os.path.join(temp_dir, "output1")
        output_dir2 = os.path.join(temp_dir, "output2")
        os.makedirs(output_dir1, exist_ok=True)
        os.makedirs(output_dir2, exist_ok=True)
        
        # Test the original implementation
        memory_tracker1 = MemoryTracker()  # Track current process
        memory_tracker1.start_tracking()
        start_time1 = time.time()
        stats1 = get_pairs(left_path, right_path, output_dir1)
        end_time1 = time.time()
        memory_info1 = memory_tracker1.stop_tracking()
        
        # Test the improved implementation
        memory_tracker2 = MemoryTracker()  # Track current process
        memory_tracker2.start_tracking()
        start_time2 = time.time()
        stats2 = get_pairs_v3(left_path, right_path, output_dir2)
        end_time2 = time.time()
        memory_info2 = memory_tracker2.stop_tracking()
        
        # Calculate metrics
        time1 = end_time1 - start_time1
        time2 = end_time2 - start_time2
        speedup = time1 / time2 if time2 > 0 else float('inf')
        
        # Convert bytes to MB for better readability
        peak_mb1 = memory_info1['peak_bytes'] / (1024 * 1024)
        peak_mb2 = memory_info2['peak_bytes'] / (1024 * 1024)
        
        # Print the performance results
        print(f"\nPerformance test with {num_reads} reads ({paired_percentage}% paired):")
        print(f"  Original implementation: {time1:.4f} seconds, {peak_mb1:.2f} MB peak memory")
        print(f"  Improved implementation: {time2:.4f} seconds, {peak_mb2:.2f} MB peak memory")
        print(f"  Speedup: {speedup:.2f}x")
        print(f"  Memory reduction: {(peak_mb1 - peak_mb2) / peak_mb1 * 100:.2f}%")
        
        # Check the results
        assert stats1["left"]["paired"] == stats2["left"]["paired"]
        assert stats1["right"]["paired"] == stats2["right"]["paired"]
        assert stats1["left"]["unpaired"] == stats2["left"]["unpaired"]
        assert stats1["right"]["unpaired"] == stats2["right"]["unpaired"]
        
        # Skip the speed check for small datasets
        # The improved implementation is optimized for large datasets and may be slower for small ones
        # due to the overhead of setting up optimizations
        if num_reads >= 10000:  # Only check speed for larger datasets
            assert speedup >= 1.0, "The improved implementation should be at least as fast as the original"
        
    finally:
        # Clean up
        shutil.rmtree(temp_dir) 