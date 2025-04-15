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

        # Positive lithiation
        Dsp_x = [-10, 0, 0.034, 0.068, 0.136, 0.204, 0.272, 0.340, 0.408,
                 0.476, 0.544, 0.612, 0.680, 1, 10]
        Dsp_y = [7e-15, 7e-15, 7e-15, 8e-15, 9e-15, 6.3e-15, 8.9e-15, 12e-15,
                 20e-15, 56e-15, 90e-15, 110e-15, 230e-15, 230e-15, 230e-15]

        self._Dsp_interp = CubicSpline(Dsp_x, Dsp_y)

        # Negative lithiation (i.e., delithiating)
        Dsn_x = [-10, 0, 0.456, 0.558, 0.660, 0.728, 0.762, 0.796, 0.830,
                 0.932, 0.976, 1, 10]
        Dsn_y = [85e-15, 85e-15, 85e-15, 45e-15, 25e-15, 12e-15, 14e-15,
                 12e-15, 10e-15, 6.0e-15, 6.0e-15, 6.0e-15, 6.0e-15]

        self._Dsn_interp = CubicSpline(Dsn_x, Dsn_y)

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

        scalar = np.isscalar(x)

        x = np.asarray(x)
        cond = np.asarray(fluxdir) >= 0

        Ds = np.empty_like(x)
        Ds[cond] = 0.0001 * self._Dsp_interp(x[cond])
        Ds[~cond] = 0.0001 * self._Dsn_interp(x[~cond])

        return Ds.item() if scalar else Ds

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
