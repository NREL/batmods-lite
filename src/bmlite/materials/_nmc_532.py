import numpy as np


class NMC532Fast:

    def __init__(self, alpha_a: float, alpha_c: float, Li_max: float) -> None:
        """
        Computationally fast NMC532 kinetic and transport properties.

        Differs from ``NMC532Slow`` because the equilibrium potential is not
        piecewise here, making it less accurate, but faster to evaluate.

        Parameters
        ----------
        alpha_a : float
            Anodic symmetry factor in Butler-Volmer expression [-].
        alpha_c : float
            Cathodic symmetry factor in Butler-Volmer expression [-].
        Li_max : float
            Maximum lithium concentration in solid phase [kmol/m^3].

        """

        self.alpha_a = alpha_a
        self.alpha_c = alpha_c
        self.Li_max = Li_max

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
            Lithium diffusivity in the solid phase [m^2/s].

        """

        from .. import Constants

        c = Constants()

        A = np.array([
            -2.509010843479270e+2,
             2.391026725259970e+3,
            -4.868420267611360e+3,
            -8.331104102921070e+1,
             1.057636028329000e+4,
            -1.268324548348120e+4,
             5.016272167775530e+3,
             9.824896659649480e+2,
            -1.502439339070900e+3,
             4.723709304247700e+2,
            -6.526092046397090e+1,
        ])

        Ds = np.exp(-30e6/c.R * (1./T - 1./303.15)) \
           * 2.25 * 10.0**(np.polyval(A, x))

        return Ds

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
            Lithium ion concentration in the local electrolyte [kmol/m^3].
        T : float
            Battery temperature [K].
        fluxdir : float | 1D array
            Lithiation direction: +1 for lithiation, -1 for delithiation, 0 for
            rest. Used for direction-dependent parameters. Ensure the zero case
            is handled explicitly or via a default (lithiating or delithiating).

        Returns
        -------
        i0 : float | 1D array
            Exchange current density [A/m^2].

        """

        from .. import Constants

        c = Constants()

        A = np.array([
             1.650452829641290e+1,
            -7.523567141488800e+1,
             1.240524690073040e+2,
            -9.416571081287610e+1,
             3.249768821737960e+1,
            -3.585290065824760e+0,
        ])

        i0 = 9.*(C_Li/1.2)**self.alpha_a * np.polyval(A, x) \
           * np.exp(-30e6/c.R * (1./T - 1./303.15))

        return i0

    def get_Eeq(self, x: float | np.ndarray) -> float | np.ndarray:
        """
        Calculate the equilibrium potential given the intercalation fraction.

        Parameters
        ----------
        x : float | 1D array
            Lithium intercalation fraction at ``r = R_s`` [-].

        Returns
        -------
        Eeq : float | 1D array
            Equilibrium potential [V].

        """

        A = np.array([
            -3.640117692001490e+3,
             1.317657544484270e+4,
            -1.455742062291360e+4,
            -1.571094264365090e+3,
             1.265630978512400e+4,
            -2.057808873526350e+3,
            -1.074374333186190e+4,
             8.698112755348720e+3,
            -8.297904604107030e+2,
            -2.073765547574810e+3,
             1.190223421193310e+3,
            -2.724851668445780e+2,
             2.723409218042130e+1,
            -4.158276603609060e+0,
             5.314735633000300e+0,
        ])

        B = np.array([
            -5.573191762723310e-4,
             6.560240842659690e+0,
             4.148209275061330e+1,
        ])

        Eeq = B[0]*np.exp(B[1] * x**B[2]) + np.polyval(A, x)

        return Eeq


class NMC532Slow(NMC532Fast):

    def __init__(self, alpha_a: float, alpha_c: float, Li_max: float) -> None:
        """
        Computationally slow NMC532 kinetic and transport properties.

        Differs from ``NMC532Fast`` because the equilibrium potential is
        piecewise here, making it more accurate, but slower to evaluate.

        Parameters
        ----------
        alpha_a : float
            Anodic symmetry factor in Butler-Volmer expression [-].
        alpha_c : float
            Cathodic symmetry factor in Butler-Volmer expression [-].
        Li_max : float
            Maximum lithium concentration in solid phase [kmol/m^3].

        """

        import os

        import pandas as pd
        from scipy.interpolate import CubicSpline

        self.alpha_a = alpha_a
        self.alpha_c = alpha_c
        self.Li_max = Li_max

        csvfile = os.path.dirname(__file__) + '/data/nmc532_ocv.csv'
        df = pd.read_csv(csvfile).sort_values(by='x')

        self.x_min = df['x'].min()
        self.x_max = df['x'].max()
        self._Eeq_spline = CubicSpline(df['x'], df['V'])

    def get_Eeq(self, x: float | np.ndarray) -> float | np.ndarray:
        """
        Calculate the equilibrium potential given the intercalation fraction.

        Parameters
        ----------
        x : float | 1D array
            Lithium intercalation fraction at ``r = R_s`` [-].

        Returns
        -------
        Eeq : float | 1D array
            Equilibrium potential [V].

        Raises
        ------
        ValueError :
            x is out of bounds [x_min, x_max].

        """

        if isinstance(x, float):
            if x < self.x_min or x > self.x_max:
                raise ValueError(f'x is out of bounds [{self.x_min},'
                                 f' {self.x_max}].')

        if not isinstance(x, float):
            if np.any(x < self.x_min) or np.any(x > self.x_max):
                raise ValueError(f'x is out of bounds [{self.x_min},'
                                 f' {self.x_max}].')

        return self._Eeq_spline(x)


class NMC532SlowExtrap(NMC532Fast):

    def __init__(self, alpha_a: float, alpha_c: float, Li_max: float) -> None:
        """
        Computationally slow NMC532 kinetic and transport properties.

        Differs from ``NMC532Slow`` because the piecewise equilibrium
        potential is extrapolated to be valid on the full [0, 1] range.

        Parameters
        ----------
        alpha_a : float
            Anodic symmetry factor in Butler-Volmer expression [-].
        alpha_c : float
            Cathodic symmetry factor in Butler-Volmer expression [-].
        Li_max : float
            Maximum lithium concentration in solid phase [kmol/m^3].

        """

        import os

        import pandas as pd
        from scipy.interpolate import CubicSpline

        self.alpha_a = alpha_a
        self.alpha_c = alpha_c
        self.Li_max = Li_max

        csvfile = os.path.dirname(__file__) + '/data/nmc532_ocv_extrap.csv'
        df = pd.read_csv(csvfile).sort_values(by='x')

        self.x_min = df['x'].min()
        self.x_max = df['x'].max()
        self._Eeq_spline = CubicSpline(df['x'], df['V'])

    def get_Eeq(self, x: float | np.ndarray) -> float | np.ndarray:
        """
        Calculate the equilibrium potential given the intercalation fraction.

        Parameters
        ----------
        x : float | 1D array
            Lithium intercalation fraction at ``r = R_s`` [-].

        Returns
        -------
        Eeq : float | 1D array
            Equilibrium potential [V].

        Raises
        ------
        ValueError :
            x is out of bounds [x_min, x_max].

        """

        return self._Eeq_spline(x)
