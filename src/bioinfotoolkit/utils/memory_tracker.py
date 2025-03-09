#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Memory Tracker

This module provides utilities for tracking memory usage of Python processes.
It consolidates functionality from multiple memory tracker modules into one.
"""

import os
import sys
import time
import platform
import logging
import threading
from typing import Dict, List, Optional, Union, Callable, Any, Tuple

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

# Platform-specific imports
if platform.system() == 'Windows':
    try:
        import psutil
    except ImportError:
        logger.warning("psutil not available. Install with: pip install psutil")
        psutil = None
elif platform.system() == 'Linux':
    try:
        import resource
    except ImportError:
        logger.warning("resource module not available.")
        resource = None
elif platform.system() == 'Darwin':  # macOS
    try:
        import resource
        import subprocess
    except ImportError:
        logger.warning("Required modules not available.")
        resource = None
        subprocess = None

class MemoryTracker:
    """
    A utility for tracking memory usage of Python processes.
    
    This class provides methods to monitor memory usage of the current process 
    or a specific function. It supports different platforms (Windows, Linux, macOS)
    and can be used as a context manager or decorator.
    """
    
    def __init__(self, 
                 interval: float = 1.0, 
                 log_level: int = logging.INFO,
                 include_children: bool = True):
        """
        Initialize the memory tracker.
        
        Args:
            interval: Polling interval in seconds
            log_level: Logging level
            include_children: Whether to include child processes in memory measurements
        """
        self.interval = interval
        self.log_level = log_level
        self.include_children = include_children
        self.measurements: List[Dict[str, Any]] = []
        self.running = False
        self.monitor_thread: Optional[threading.Thread] = None
        self.start_time = 0.0
        self.end_time = 0.0
        self.peak_memory = 0
        self.current_memory = 0
        
        # Set up platform-specific measuring method
        self.measure_memory = self._get_platform_measure_method()
    
    def _get_platform_measure_method(self) -> Callable[[], int]:
        """
        Get the appropriate memory measurement method for the current platform.
        
        Returns:
            A function that returns memory usage in bytes
        """
        system = platform.system()
        
        if system == 'Windows':
            if psutil:
                return self._measure_memory_windows
            else:
                logger.warning("Using fallback memory measurement on Windows. Install psutil for accurate tracking.")
                return self._measure_memory_fallback
        
        elif system == 'Linux':
            if resource:
                return self._measure_memory_linux
            else:
                return self._measure_memory_fallback
        
        elif system == 'Darwin':  # macOS
            if resource and subprocess:
                return self._measure_memory_macos
            else:
                return self._measure_memory_fallback
        
        else:
            logger.warning(f"Unsupported platform: {system}. Using fallback memory measurement.")
            return self._measure_memory_fallback
    
    def _measure_memory_windows(self) -> int:
        """Measure memory usage on Windows."""
        process = psutil.Process(os.getpid())
        memory_info = process.memory_info()
        
        if self.include_children:
            # Include memory of child processes
            mem_total = memory_info.rss
            for child in process.children(recursive=True):
                try:
                    mem_total += child.memory_info().rss
                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    pass
            return mem_total
        else:
            return memory_info.rss
    
    def _measure_memory_linux(self) -> int:
        """Measure memory usage on Linux."""
        # Get resident set size (RSS) from /proc/self/statm
        with open('/proc/self/statm', 'r') as f:
            fields = f.read().split()
            rss = int(fields[1]) * os.sysconf('SC_PAGE_SIZE')
        
        # Include child processes if requested
        if self.include_children and psutil:
            process = psutil.Process(os.getpid())
            for child in process.children(recursive=True):
                try:
                    rss += child.memory_info().rss
                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    pass
        
        return rss
    
    def _measure_memory_macos(self) -> int:
        """Measure memory usage on macOS."""
        # Use resource module for current process
        rusage = resource.getrusage(resource.RUSAGE_SELF)
        memory = rusage.ru_maxrss * 1024  # macOS reports in KB
        
        # Include memory of child processes if requested
        if self.include_children:
            rusage_children = resource.getrusage(resource.RUSAGE_CHILDREN)
            memory += rusage_children.ru_maxrss * 1024
        
        return memory
    
    def _measure_memory_fallback(self) -> int:
        """Fallback memory measurement method when platform-specific methods aren't available."""
        try:
            import tracemalloc
            if not tracemalloc.is_tracing():
                tracemalloc.start()
            
            snapshot = tracemalloc.take_snapshot()
            return sum(trace.size for trace in snapshot.traces)
        except ImportError:
            logger.warning("tracemalloc not available. Memory tracking will be inaccurate.")
            return 0
    
    def start(self) -> None:
        """Start memory tracking."""
        if self.running:
            logger.warning("Memory tracker is already running.")
            return
        
        self.running = True
        self.measurements = []
        self.start_time = time.time()
        self.peak_memory = 0
        
        # Start the monitoring thread
        self.monitor_thread = threading.Thread(
            target=self._monitor_memory,
            daemon=True
        )
        self.monitor_thread.start()
    
    def stop(self) -> Tuple[float, int]:
        """
        Stop memory tracking.
        
        Returns:
            Tuple of (elapsed_time, peak_memory)
        """
        if not self.running:
            logger.warning("Memory tracker is not running.")
            return 0.0, 0
        
        self.running = False
        if self.monitor_thread:
            self.monitor_thread.join(timeout=self.interval*2)
        
        self.end_time = time.time()
        elapsed_time = self.end_time - self.start_time
        
        return elapsed_time, self.peak_memory
    
    def _monitor_memory(self) -> None:
        """Monitor memory usage in a background thread."""
        while self.running:
            try:
                memory = self.measure_memory()
                self.current_memory = memory
                self.peak_memory = max(self.peak_memory, memory)
                
                timestamp = time.time() - self.start_time
                self.measurements.append({
                    'timestamp': timestamp,
                    'memory': memory
                })
                
                logger.log(self.log_level, f"Memory usage: {self._format_bytes(memory)}")
                
            except Exception as e:
                logger.error(f"Error measuring memory: {e}")
            
            # Wait for the next interval
            time.sleep(self.interval)
    
    def _format_bytes(self, bytes_value: int) -> str:
        """
        Format bytes as a human-readable string.
        
        Args:
            bytes_value: Number of bytes
            
        Returns:
            Formatted string (e.g. "4.5 MB")
        """
        units = ['B', 'KB', 'MB', 'GB', 'TB']
        size = float(bytes_value)
        unit_index = 0
        
        while size >= 1024.0 and unit_index < len(units) - 1:
            size /= 1024.0
            unit_index += 1
        
        return f"{size:.2f} {units[unit_index]}"
    
    def get_summary(self) -> Dict[str, Any]:
        """
        Get a summary of memory usage.
        
        Returns:
            Dictionary with memory usage statistics
        """
        elapsed_time = self.end_time - self.start_time if self.end_time else time.time() - self.start_time
        
        return {
            'peak_memory': self.peak_memory,
            'peak_memory_formatted': self._format_bytes(self.peak_memory),
            'elapsed_time': elapsed_time,
            'elapsed_time_formatted': f"{elapsed_time:.2f} seconds",
            'measurements': self.measurements
        }
    
    def __enter__(self) -> 'MemoryTracker':
        """Start tracking when used as a context manager."""
        self.start()
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        """Stop tracking when exiting context manager."""
        self.stop()
        
        # Log peak memory usage
        logger.log(self.log_level, f"Peak memory usage: {self._format_bytes(self.peak_memory)}")
        logger.log(self.log_level, f"Elapsed time: {self.end_time - self.start_time:.2f} seconds")
    
    def __call__(self, func: Callable) -> Callable:
        """Use as a decorator to track memory usage of a function."""
        def wrapper(*args, **kwargs):
            with self:
                result = func(*args, **kwargs)
            return result
        
        return wrapper

# Simple functions for quick memory monitoring
def get_current_memory() -> int:
    """
    Get current memory usage of the process.
    
    Returns:
        Current memory usage in bytes
    """
    tracker = MemoryTracker()
    return tracker.measure_memory()

def get_current_memory_formatted() -> str:
    """
    Get current memory usage as a formatted string.
    
    Returns:
        Formatted string (e.g. "4.5 MB")
    """
    tracker = MemoryTracker()
    memory = tracker.measure_memory()
    return tracker._format_bytes(memory)

def track_memory(func=None, *, interval: float = 1.0, log_level: int = logging.INFO):
    """
    Decorator to track memory usage of a function.
    
    Args:
        func: Function to decorate
        interval: Polling interval in seconds
        log_level: Logging level
        
    Returns:
        Decorated function
    """
    if func is None:
        return lambda f: MemoryTracker(interval=interval, log_level=log_level)(f)
    return MemoryTracker(interval=interval, log_level=log_level)(func)

# Export public API
__all__ = [
    'MemoryTracker',
    'measure_memory_usage',
    'measure_script_memory',
    'format_bytes'
] 