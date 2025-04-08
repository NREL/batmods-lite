"""
Material Properties Package
---------------------------
The ``materials`` package contains kinetic and transport properties for common
battery materials, including both electrode active materials and electrolytes.

"""

from ._gen2_electrolyte import Gen2Electrolyte

from ._graphite import GraphiteFast, GraphiteSlow, GraphiteSlowExtrap

from ._lfp import LFPInterp
from ._nmc_532 import NMC532Fast, NMC532Slow, NMC532SlowExtrap

__all__ = [
    Gen2Electrolyte,
    GraphiteFast,
    GraphiteSlow,
    GraphiteSlowExtrap,
    LFPInterp,
    NMC532Fast,
    NMC532Slow,
    NMC532SlowExtrap,
]
