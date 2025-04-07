from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:  # pragma: no cover
    from numpy import ndarray

    from .domains import Electrode


class Hysteresis:

    def __init__(self, domain: Electrode, **options) -> None:
        self.domain = domain
        domain.g_hyst = options.pop('g_hyst')
        domain.M_hyst = options.pop('M_hyst')
        domain.hyst0 = options.pop('hyst0')

    def make_mesh(self, pshift: int = 0) -> None:

        domain = self.domain
        if domain._name == 'anode':
            domain.ptr['hyst'] = domain.ptr['phi_ed'] + 1

        elif domain._name == 'cathode':
            for k in domain.ptr.keys():
                domain.ptr[k] += 1
            domain.ptr['hyst'] = 0 + pshift

    def sv0(self, sv0: ndarray) -> None:

        domain = self.domain
        start = domain.ptr['start']
        sv0[domain.ptr['hyst'] - start] = domain.hyst0

    def algidx(self, algdix: ndarray) -> None:
        pass

    def to_dict(self, sol: object) -> None:
        domain = self.domain
        hyst = sol.y[:, domain.ptr['hyst']]
        return {'hyst': hyst}
