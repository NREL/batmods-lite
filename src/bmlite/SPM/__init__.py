"""
Single Particle Model Package
-----------------------------
A packaged single particle model (SPM). Build a model using the ``Simulation``
class, and run any available experiments using its "run" methods, e.g.
``run_CC()``. The experiments return ``Solution`` class instances with post
processing, plotting, and saving methods.

"""

from ._simulation import Simulation
from ._solutions import StepSolution, CycleSolution

from . import dae
from . import domains
from . import postutils
from . import submodels

__all__ = [
    'Simulation',
    'StepSolution',
    'CycleSolution',
    'dae',
    'domains',
    'postutils',
    'submodels',
]
