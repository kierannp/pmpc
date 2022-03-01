"""
Unit and regression test for the pmpc package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import pmpc


def test_pmpc_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "pmpc" in sys.modules
