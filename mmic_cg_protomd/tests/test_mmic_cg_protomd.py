"""
Unit and regression test for the mmic_cg_protomd package.
"""

# Import package, test suite, and other packages as needed
import mmelemental
from mmelemental.models import Molecule
from mmic_cg.models import CoarseInput, CoarseOutput
from mmic_cg_protomd.components import CoarseProtoMDComponent

import pytest
import sys
import os
import pytest
import mm_data


def test_mmic_cg_protomd_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "mmic_cg_protomd" in sys.modules


def test_compute_component():
    """
    This test is used to check
    if protomd can be successfully used
    to coarse=grain a system. The method here
    is SpaceWarping and a water molecule is used
    """

    mol = mmelemental.models.Molecule.from_file(mm_data.mols["water-mol.json"])

    inputs = CoarseInput(
        molecule={"mol": mol},
        method="spacewarping",
        schema_name="test",
        schema_version=1.0,
        keywords={"kmax": 1},
    )

    outputs = CoarseProtoMDComponent.compute(inputs)
