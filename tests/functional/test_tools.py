"""Functional tests for the tools scripts."""

import os
import pytest
import subprocess
import sys
import tempfile
from pathlib import Path


def test_script_template():
    """Test the script_template script."""
    # Get script path
    script_path = Path(__file__).parent.parent.parent / "src" / "scripts" / "tools" / "script_template.py"
    
    # Run the script with help option to verify it works
    cmd = [
        sys.executable,
        str(script_path),
        "-h"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    # Check that the script executed successfully
    assert result.returncode == 0
    # Just verify the help output contains some expected text
    assert "usage:" in result.stdout

