"""
BATMODS-lite
============
Battery Analysis and Training Models for Optimization and Degradation Studies
(BATMODS) is a Python package with an API for pre-built battery models. The
original purpose of the package was to quickly generate synthetic data for
machine learning models to train with. However, the models are generally useful
for any battery simulation or analysis. BATMODS-lite includes the following:

1) A library and API for pre-built battery models
2) Kinetic/transport properties for common battery materials

How to use the documentation
----------------------------
Documentation is accessible via Python's ``help()`` function which prints
docstrings from the specified package, module, function, class, etc. (e.g.,
``help(bmlite.SPM)``). In addition, you can access the documentation
by calling the built-in ``bmlite.docs()`` method to open a locally run website.
The website includes search functionality and examples, beyond the code
docstrings.

Viewing documentation using IPython
-----------------------------------
Start IPython and import ``bmlite`` under the alias ``bm``. To see what's
available in ``bmlite``, type ``bm.<TAB>`` (where ``<TAB>`` refers to the TAB
key). To view the type hints and brief descriptions, type an open parenthesis
``(`` after any function, method, class, etc. (e.g., ``bm.Constants(``).

"""

# Core package
from ._core import (
    Constants,
    Experiment,
    IDASolver,
    IDAResult,
    _templates,
)

# Model subpackages
from . import P2D
from . import SPM

# Other submodules/subpackages
from . import mathutils
from . import materials
from . import mesh
from . import plotutils

__version__ = '0.0.2'

__all__ = [
    'Constants',
    'Experiment',
    'IDASolver',
    'IDAResult',
    '_templates',
    'P2D',
    'SPM',
    'mathutils',
    'materials',
    'mesh',
    'plotutils',
]
