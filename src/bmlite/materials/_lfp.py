from numbers import Real

import numpy as np


class LFPInterp:

    def __init__(self, alpha_a: float, alpha_c: float, Li_max: float) -> None:
        """
        LFP kinetic and transport properties, interpolated from data.

        Parameters
        ----------
        alpha_a : float
            Anodic symmetry factor in Butler-Volmer expression [-].
        alpha_c : float
            Cathodic symmetry factor in Butler-Volmer expression [-].
        Li_max : float
            Maximum lithium concentration in solid phase [kmol/m3].

        """

        import os

        import pandas as pd
        from scipy.interpolate import CubicSpline

        self.alpha_a = alpha_a
        self.alpha_c = alpha_c
        self.Li_max = Li_max

        # Diffusivity polynomial coeffs, for Ds = 10**P(x)
        self._Ds_coeffs = np.array([
             349.29800720, -925.34470810, 875.56218099, -268.24357721,
            -102.31082863,  100.13715321, -27.74353094, -17.700056810,
        ])

        # OCV and hysteresis, from csv data file
        csvfile = os.path.dirname(__file__) + '/data/lfp_ocv.csv'
        df = pd.read_csv(csvfile).sort_values(by='x')

        self._Eeq_interp = CubicSpline(df['x'], df['V_avg'])
        self._Mhyst_interp = CubicSpline(df['x'], df['M_hyst'])

    def get_Ds(self, x: float | np.ndarray, T: float,
               fluxdir: float | np.ndarray) -> float | np.ndarray:
        """
        Calculate the lithium diffusivity in the solid phase given the local
        intercalation fraction ``x`` and temperature ``T``.

        Parameters
        ----------
        x : float | 1D array
            Lithium intercalation fraction [-].
        T : float
            Battery temperature [K].
        fluxdir : float | 1D array
            Lithiation direction: +1 for lithiation, -1 for delithiation, 0 for
            rest. Used for direction-dependent parameters. Ensure the zero case
            is handled explicitly or via a default (lithiating or delithiating).

        Returns
        -------
        Ds : float | 1D array
            Lithium diffusivity in the solid phase [m2/s].

        """

        Ds = 4.014e-17*np.ones_like(x)  # 10**np.polyval(self._Ds_coeffs, x)

        return Ds.item() if np.isscalar(x) else Ds

    def get_i0(self, x: float | np.ndarray, C_Li: float | np.ndarray,
               T: float, fluxdir: float | np.ndarray) -> float | np.ndarray:
        """
        Calculate the exchange current density given the intercalation
        fraction ``x`` at the particle surface, the local lithium ion
        concentration ``C_Li``, and temperature ``T``. The input types for
        ``x`` and ``C_Li`` should both be the same (i.e., both float or both
        1D arrays).

        Parameters
        ----------
        x : float | 1D array
            Lithium intercalation fraction at ``r = R_s`` [-].
        C_Li : float | 1D array
            Lithium ion concentration in the local electrolyte [kmol/m3].
        T : float
            Battery temperature [K].
        fluxdir : float | 1D array
            Lithiation direction: +1 for lithiation, -1 for delithiation, 0 for
            rest. Used for direction-dependent parameters. Ensure the zero case
            is handled explicitly or via a default (lithiating or delithiating).

        Returns
        -------
        i0 : float | 1D array
            Exchange current density [A/m2].

        """

        from .. import Constants

        c = Constants()

        # Avoid floating point errors
        if isinstance(x, Real):
            if (x < 0 and self.alpha_c < 1) or (x > 1 and self.alpha_a < 1):
                raise ValueError('x is out of [0, 1] during i0 calculation')
        elif isinstance(x, np.ndarray):
            if ((any(x.flatten() < 0) and self.alpha_c < 1)
               or (any(x.flatten() > 1) and self.alpha_a < 1)):
                raise ValueError('x is out of [0, 1] during i0 calculation')

        i0 = 0.27 * np.exp(-30e6 / c.R * (1 / T - 1 / 303.15)) \
           * C_Li**self.alpha_a * (self.Li_max * x)**self.alpha_c \
           * (self.Li_max - self.Li_max * x)**self.alpha_a

        return i0

    def get_Eeq(self, x: float | np.ndarray) -> float | np.ndarray:
        """
        Calculate the equilibrium potential given the surface intercalation
        fraction ``x`` at the particle surface.

        Parameters
        ----------
        x : float
            Lithium intercalation fraction at ``r = R_s`` [-].

        Returns
        -------
        Eeq : float
            Equilibrium potential [V].

        """

        return self._Eeq_interp(x)

    def get_Mhyst(self, x: float | np.ndarray) -> float | np.ndarray:
        """
        Calculate the hysteresis magnitude given the surface intercalation
        fraction ``x`` at the particle surface.

        Parameters
        ----------
        x : float | 1D array
            Lithium intercalation fraction at ``r = R_s`` [-].

        Returns
        -------
        M_hyst : float | 1D array
            Hysteresis magnitude [V].

        """

        return self._Mhyst_interp(x)
