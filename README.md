# Coarse-Graining Component

[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/RlyehAD/mmic_cg_protomd/workflows/CI/badge.svg)](https://github.com/RlyehAD/mmic_cg_protomd/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/RlyehAD/mmic_cg_protoMD/branch/master/graph/badge.svg)](https://codecov.io/gh/RlyehAD/mmic_cg_protoMD/branch/master)
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/RlyehAD/mmic_cg_protomd.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/RlyehAD/mmic_cg_protomd/context:python)


This is part of the [MolSSI](http://molssi.org) Molecular Mechanics Interoperable Components ([MMIC](https://github.com/MolSSI/mmic)) project. This package provides a component for a tactic component for [mmic_cg](https://github.com/MolSSI/mmic_cg.git) using the [ProtoMD](https://github.com/CTCNano/proto_md.git) software suite.

# Basic Usage

## Preparing Input
```python
# Import main component for coarse-graining
from mmic_cg_protomd.components import CoarseProtoMDComponent

# Import a molecule model that complies with MMSchema
from mmelemental.models import Molecule
from mmic_cg.models import InputCoarse

# Create MMSchema molecule
mol = Molecule.from_file(path_to_file)

# Create input for coarse-graining a molecule with protoMD
cgInp = InputCoarse(
    "molecule"=mol, 
    "method"="spacewarping",
    "schema_name"="test",
    "schema_version"=1.0,
    "method_keywords"={
        "kmax": 0,
    },
)
```

## Run The CG Algorithm
```python
# Execute coarse-graining
cgOut = CoarseProtoMDComponent.compute(cgInp)
```

# Extract MMSchema CG mol
```python
cgMol = cgOut.mol
cgVel = cgOut.mol["velocities"]
cgPos = cgOut.mol["geometry"]
```

### Copyright

Copyright (c) 2021, Xu Guo, Andrew Abi-Mansour


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.1.
