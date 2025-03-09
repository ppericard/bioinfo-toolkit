#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Bioinfo Logger

Description: Smart logging utility for bioinformatics tools that provides
             adaptive reporting based on dataset size.
"""

import logging
import time
import sys
from typing import Dict, Any, Optional, Callable

# Configure default logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)


class BioinfLogger:
    """
    A logger that adapts its verbosity based on dataset size.
    Provides smart progress reporting and performance metrics.
    """
    def __init__(self, name: str, level: int = logging.INFO):
        """
        Initialize the logger.
        
        Args:
            name: Logger name (typically the script name)
            level: Initial logging level
        """
        self.logger = logging.getLogger(name)
        self.logger.setLevel(level)
        
        # Performance tracking attributes
        self.start_time = time.time()
        self.last_progress_time = self.start_time
        self.progress_interval = 0.5  # seconds between progress updates
        
        # Adaptive reporting attributes
        self.record_count = 0
        self.last_report_count = 0
        self.report_frequency = 100000  # Initial reporting frequency

    def set_level(self, level: int) -> None:
        """Set the logging level."""
        self.logger.setLevel(level)
    
    def debug(self, message: str) -> None:
        """Log a debug message."""
        self.logger.debug(message)
    
    def info(self, message: str) -> None:
        """Log an info message."""
        self.logger.info(message)
    
    def warning(self, message: str) -> None:
        """Log a warning message."""
        self.logger.warning(message)
    
    def error(self, message: str) -> None:
        """Log an error message."""
        self.logger.error(message)
    
    def critical(self, message: str) -> None:
        """Log a critical message."""
        self.logger.critical(message)
    
    def start_progress(self, message: str) -> None:
        """
        Start a progress tracking operation.
        
        Args:
            message: Initial progress message
        """
        self.start_time = time.time()
        self.last_progress_time = self.start_time
        self.record_count = 0
        self.last_report_count = 0
        self.info(f"Starting: {message}")
    
    def update_progress(self, count: int, stats: Optional[Dict[str, Any]] = None, 
                       force: bool = False) -> None:
        """
        Update progress for an operation, with adaptive reporting frequency.
        
        Args:
            count: Current record count
            stats: Optional statistics to include in the report
            force: If True, always report regardless of timing
        """
        now = time.time()
        self.record_count = count
        
        # Determine if we should report progress
        should_report = force or (
            # Enough time has passed
            (now - self.last_progress_time >= self.progress_interval) and 
            # Enough new records processed
            (self.record_count - self.last_report_count >= self.report_frequency)
        )
        
        if should_report:
            elapsed = now - self.start_time
            rate = self.record_count / elapsed if elapsed > 0 else 0
            
            # Format count with commas
            formatted_count = f"{self.record_count:,}"
            
            # Basic progress message
            message = f"Processed {formatted_count} records ({rate:.2f}/sec)"
            
            # Add statistics if provided
            if stats:
                stat_parts = []
                for key, value in stats.items():
                    if isinstance(value, (int, float)):
                        stat_parts.append(f"{key}: {value:,}")
                    else:
                        stat_parts.append(f"{key}: {value}")
                
                if stat_parts:
                    message += f" - {', '.join(stat_parts)}"
            
            self.logger.info(message)
            
            # Update tracking variables
            self.last_progress_time = now
            self.last_report_count = self.record_count
            
            # Adapt reporting frequency based on processing rate
            # For faster processing, report less frequently
            if rate > 100000:  # Very fast processing
                self.report_frequency = 1000000
                self.progress_interval = 5.0
            elif rate > 10000:  # Fast processing
                self.report_frequency = 100000
                self.progress_interval = 2.0
            elif rate > 1000:  # Medium processing
                self.report_frequency = 10000
                self.progress_interval = 1.0
            else:  # Slow processing
                self.report_frequency = 1000
                self.progress_interval = 0.5
    
    def finish_progress(self, message: str = "Completed") -> Dict[str, float]:
        """
        Finish progress tracking and report final statistics.
        
        Args:
            message: Message to display on completion
            
        Returns:
            Dictionary with performance statistics
        """
        end_time = time.time()
        elapsed = end_time - self.start_time
        rate = self.record_count / elapsed if elapsed > 0 else 0
        
        stats = {
            "elapsed_seconds": elapsed,
            "records_processed": self.record_count,
            "records_per_second": rate
        }
        
        self.logger.info(f"{message} in {elapsed:.2f} seconds ({rate:.2f} records/second)")
        return stats
    
    def get_progress_callback(self, message_prefix: str = "") -> Callable:
        """
        Return a callback function that updates progress.
        
        Args:
            message_prefix: Optional prefix for progress messages
            
        Returns:
            Callback function for progress updates
        """
        prefix = f"{message_prefix}: " if message_prefix else ""
        
        def callback(count, stats=None):
            if count > self.record_count:
                self.update_progress(count, stats)
        
        return callback


def get_logger(name: str, level: int = logging.INFO) -> BioinfLogger:
    """
    Factory function to create a BioinfLogger.
    
    Args:
        name: Logger name
        level: Initial logging level
        
    Returns:
        Configured BioinfLogger instance
    """
    return BioinfLogger(name, level) 