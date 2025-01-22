"""
Pseudo-2D Model Package
-----------------------
A packaged pseudo-2D (P2D) model. Build a model using the ``Simulation`` class,
and run any available experiments using its "run" methods, e.g. ``run_CC()``.
The experiments return ``Solution`` class instances with post processing,
plotting, and saving methods.

"""

from ._simulation import Simulation
from ._solutions import StepSolution, CycleSolution

from . import dae
from . import domains
from . import postutils

__all__ = [
    'Simulation',
    'StepSolution',
    'CycleSolution',
    'dae',
    'domains',
    'postutils',
]
