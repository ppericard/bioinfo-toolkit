"""Utility modules for bioinfo-toolkit."""

from bioinfotoolkit.utils.bioinfo_logger import get_logger, BioinfLogger
from bioinfotoolkit.utils.fastq_utils import (
    open_file,
    extract_read_id,
    read_fastq_records,
    chunk_iterator,
    build_read_id_set,
    process_fastq_file
)
from bioinfotoolkit.utils.memory_tracker import (
    MemoryTracker,
    get_current_memory,
    get_current_memory_formatted,
    track_memory
)

__all__ = [
    # Logging
    'get_logger', 'BioinfLogger',
    
    # FASTQ utilities
    'open_file', 'extract_read_id', 'read_fastq_records',
    'chunk_iterator', 'build_read_id_set', 'process_fastq_file',
    
    # Memory tracking
    'MemoryTracker', 'get_current_memory', 
    'get_current_memory_formatted', 'track_memory'
] 