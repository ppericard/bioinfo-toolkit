#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Memory Tracker

Description: Utilities for tracking and reporting memory usage in Python scripts
"""

import os
import sys
import subprocess
import tempfile
import platform
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union, Callable

# Try importing platform-specific modules
try:
    import psutil
    HAVE_PSUTIL = True
except ImportError:
    HAVE_PSUTIL = False

try:
    import resource
    HAVE_RESOURCE = True
except ImportError:
    HAVE_RESOURCE = False


class MemoryTracker:
    """
    Track memory usage during script execution using various methods.
    Supports multiple tracking mechanisms for cross-platform compatibility.
    """
    
    def __init__(self, pid: Optional[int] = None):
        """
        Initialize the memory tracker.
        
        Args:
            pid: Process ID to track, or None to track the current process
        """
        self.pid = pid or os.getpid()
        self.platform = platform.system()
        self.tracking_methods = self._get_available_tracking_methods()
        self.peak_memory = 0  # in bytes
        self.start_memory = self._get_current_memory()
        self.timeline = []  # [(timestamp, memory_usage)]
        self.is_tracking = False
        
    def _get_available_tracking_methods(self) -> List[str]:
        """
        Determine which memory tracking methods are available on this system.
        
        Returns:
            List of tracking method names
        """
        methods = []
        
        # psutil is cross-platform
        if HAVE_PSUTIL:
            methods.append('psutil')
        
        # resource module for Unix-like systems
        if HAVE_RESOURCE and self.platform != 'Windows':
            methods.append('resource')
        
        # Operating system specific methods
        if self.platform == 'Linux':
            methods.append('proc')
        elif self.platform == 'Windows':
            # Windows Performance Counters
            methods.append('wmic')
        
        return methods
    
    def _get_memory_psutil(self) -> int:
        """Get memory usage via psutil in bytes."""
        if not HAVE_PSUTIL:
            return 0
        
        try:
            process = psutil.Process(self.pid)
            # Return RSS (Resident Set Size) in bytes
            return process.memory_info().rss
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            return 0
    
    def _get_memory_resource(self) -> int:
        """Get memory usage via resource module in bytes."""
        if not HAVE_RESOURCE:
            return 0
        
        try:
            # Return peak memory (maximum resident set size) in bytes
            return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024
        except Exception:
            return 0
    
    def _get_memory_proc(self) -> int:
        """Get memory from /proc/[pid]/status on Linux systems."""
        if self.platform != 'Linux':
            return 0
        
        try:
            with open(f'/proc/{self.pid}/status', 'r') as f:
                for line in f:
                    if line.startswith('VmRSS:'):
                        # Extract value in kB and convert to bytes
                        return int(line.split()[1]) * 1024
        except (FileNotFoundError, IOError, ValueError):
            return 0
        
        return 0
    
    def _get_memory_wmic(self) -> int:
        """Get memory usage on Windows via WMIC."""
        if self.platform != 'Windows':
            return 0
        
        try:
            cmd = f'wmic process where ProcessId={self.pid} get WorkingSetSize'
            output = subprocess.check_output(cmd, shell=True, text=True)
            lines = output.strip().split('\n')
            if len(lines) >= 2:
                # Parse the second line (first line is the header)
                return int(lines[1].strip())
        except (subprocess.SubprocessError, ValueError, IndexError):
            return 0
        
        return 0
    
    def _get_current_memory(self) -> int:
        """
        Get current memory usage using the best available method.
        
        Returns:
            Memory usage in bytes
        """
        memory_usage = 0
        
        # Try each method in order of accuracy
        if 'psutil' in self.tracking_methods:
            memory_usage = max(memory_usage, self._get_memory_psutil())
        
        if 'resource' in self.tracking_methods:
            memory_usage = max(memory_usage, self._get_memory_resource())
        
        if 'proc' in self.tracking_methods:
            memory_usage = max(memory_usage, self._get_memory_proc())
        
        if 'wmic' in self.tracking_methods:
            memory_usage = max(memory_usage, self._get_memory_wmic())
        
        return memory_usage
    
    def start_tracking(self, interval: float = 0.1) -> None:
        """
        Start tracking memory usage in a background thread.
        
        Args:
            interval: Time between memory measurements in seconds
        """
        if HAVE_PSUTIL:
            self.is_tracking = True
            # Reset statistics
            self.peak_memory = self._get_current_memory()
            self.timeline = [(time.time(), self.peak_memory)]
            
            # Start monitoring in a background thread
            import threading
            
            def _monitor():
                while self.is_tracking:
                    current_memory = self._get_current_memory()
                    self.peak_memory = max(self.peak_memory, current_memory)
                    self.timeline.append((time.time(), current_memory))
                    time.sleep(interval)
            
            self.monitor_thread = threading.Thread(target=_monitor, daemon=True)
            self.monitor_thread.start()
        else:
            # If psutil isn't available, just record the start memory
            self.peak_memory = self._get_current_memory()
            self.timeline = [(time.time(), self.peak_memory)]
    
    def stop_tracking(self) -> Dict[str, Union[int, float]]:
        """
        Stop tracking memory usage and return statistics.
        
        Returns:
            Dictionary with memory usage statistics
        """
        # Measure final memory usage
        current = self._get_current_memory()
        self.peak_memory = max(self.peak_memory, current)
        self.timeline.append((time.time(), current))
        
        # Stop the monitoring thread if it exists
        self.is_tracking = False
        if hasattr(self, 'monitor_thread') and self.monitor_thread.is_alive():
            self.monitor_thread.join(timeout=1.0)
        
        # Calculate statistics
        initial = self.timeline[0][1] if self.timeline else 0
        final = self.timeline[-1][1] if self.timeline else current
        
        return {
            'peak_bytes': self.peak_memory,
            'initial_bytes': initial,
            'final_bytes': final,
            'difference_bytes': final - initial
        }


def measure_memory_usage(func: Callable, *args, **kwargs) -> Tuple[object, Dict[str, Union[int, float]]]:
    """
    Measure the memory usage of a function.
    
    Args:
        func: The function to measure
        *args: Arguments to pass to the function
        **kwargs: Keyword arguments to pass to the function
        
    Returns:
        Tuple of (function result, memory statistics)
    """
    tracker = MemoryTracker()
    tracker.start_tracking()
    
    try:
        result = func(*args, **kwargs)
    finally:
        stats = tracker.stop_tracking()
    
    return result, stats


def measure_script_memory(cmd: List[str]) -> Tuple[str, Dict[str, Union[int, float]]]:
    """
    Measure the memory usage of an external script.
    
    Args:
        cmd: Command to execute as a list of strings
        
    Returns:
        Tuple of (command output, memory statistics)
    """
    # Use built-in tools if available
    if platform.system() == 'Linux':
        # Try using /usr/bin/time -v
        try:
            # Create a temporary file for output
            with tempfile.NamedTemporaryFile(mode='w+', delete=False) as f:
                time_file = f.name
            
            # Create a temporary file for stderr
            with tempfile.NamedTemporaryFile(mode='w+', delete=False) as f:
                stderr_file = f.name
            
            # Add time command to measure memory
            time_cmd = ['/usr/bin/time', '-v', '-o', time_file]
            
            # Execute the command
            full_cmd = time_cmd + cmd
            output = subprocess.check_output(full_cmd, stderr=open(stderr_file, 'w'), text=True)
            
            # Parse time output
            memory_kb = 0
            with open(time_file, 'r') as f:
                for line in f:
                    if 'Maximum resident set size' in line:
                        # Value is in KB
                        try:
                            memory_kb = int(line.split(':')[1].strip())
                            break
                        except (ValueError, IndexError):
                            pass
            
            # Clean up temporary files
            try:
                os.unlink(time_file)
                os.unlink(stderr_file)
            except OSError:
                pass
            
            return output, {'peak_bytes': memory_kb * 1024}
            
        except (subprocess.SubprocessError, FileNotFoundError):
            # Fall back to psutil method
            pass
    
    # Default method using psutil
    if HAVE_PSUTIL:
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        
        # Track memory in the child process
        tracker = MemoryTracker(process.pid)
        tracker.start_tracking()
        
        # Wait for process to complete
        stdout, stderr = process.communicate()
        
        # Stop tracking and get stats
        stats = tracker.stop_tracking()
        
        return stdout, stats
    
    # If no method is available, just run the command without tracking
    output = subprocess.check_output(cmd, text=True)
    return output, {'peak_bytes': 0}


def format_bytes(bytes_value: int) -> str:
    """
    Format bytes into a human-readable string.
    
    Args:
        bytes_value: Value in bytes
        
    Returns:
        Formatted string (e.g., "1.23 MB")
    """
    if bytes_value < 1024:
        return f"{bytes_value} B"
    elif bytes_value < 1024 * 1024:
        return f"{bytes_value / 1024:.2f} KB"
    elif bytes_value < 1024 * 1024 * 1024:
        return f"{bytes_value / (1024 * 1024):.2f} MB"
    else:
        return f"{bytes_value / (1024 * 1024 * 1024):.2f} GB" 