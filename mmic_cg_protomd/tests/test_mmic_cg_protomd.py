"""
Unit and regression test for the mmic_cg_protomd package.
"""

# Import package, test suite, and other packages as needed
import mmic_cg_protomd
import pytest
import sys


def test_mmic_cg_protomd_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "mmic_cg_protomd" in sys.modules
