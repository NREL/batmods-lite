from numbers import Real

import numpy as np


class GraphiteFast:

    def __init__(self, alpha_a: float, alpha_c: float, Li_max: float) -> None:
        """
        Computationally fast Graphite kinetic and transport properties.

        Differs from ``GraphiteSlow`` because the equilibrium potential is
        not piecewise here, making it less accurate, but faster to evaluate.

        Parameters
        ----------
        alpha_a : float
            Anodic symmetry factor in Butler-Volmer expression [-].
        alpha_c : float
            Cathodic symmetry factor in Butler-Volmer expression [-].
        Li_max : float
            Maximum lithium concentration in solid phase [kmol/m3].

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
            Lithium diffusivity in the solid phase [m2/s].

        """

        from .. import Constants

        c = Constants()

        Ds = 3e-14 * np.exp(-30e6 / c.R * (1 / T - 1 / 303.15))

        if np.atleast_1d(x).size > 1:
            Ds = Ds * np.ones_like(x)

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

        i0 = 2.5 * 0.27 * np.exp(-30e6 / c.R * (1 / T - 1 / 303.15)) \
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

        x = np.atleast_1d(x)

        A = np.array([
        -1.059423355572770e-2, -1.453708425609560e-2, 9.089868397988610e-5,
         2.443615203087110e-2, -5.464261369950400e-1, 6.270508166379020e-1,
        -1.637520788053810e-2, -5.639025014475490e-1, 7.053886409518520e-2,
        -6.542365622896410e-2, -5.960370524233590e-1, 1.409966536648620e+0,
        -4.173226059293490e-2, -1.787670587868640e-1, 7.693844911793470e-2,
        -4.792178163846890e-1,  3.845707852011820e-3, 4.112633446959460e-2,
         6.594735004847470e-1,
        ])

        B = np.array([
        -4.364293924074990e-2, -9.449231893318330e-2, -2.046776012570780e-2,
        -8.241166396760410e-2, -7.746685789572230e-2,  3.593817905677970e-2,
        ])

        C = np.array([
        -1.731504647676420e+2,  8.252008712749000e+1,  1.233160814852810e+2,
         5.913206621637760e+1,  3.322960033709470e+1,  3.437968012320620e+0,
        -6.906367679257650e+1, -1.228217254296760e+1, -5.037944982759270e+1,
        ])

        D = np.array([
         1.059423355572770e-2, -1.453708425609560e-2, 9.089868397988610e-5,
         2.443615203087110e-2, -5.464261369950400e-1, 6.270508166379020e-1,
        -1.637520788053810e-2, -5.639025014475490e-1, 7.053886409518520e-2,
        -6.542365622896410e-2, -5.960370524233590e-1, 1.409966536648620e+0,
        -4.173226059293490e-2, -1.787670587868640e-1, 7.693844911793470e-2,
        -4.792178163846890e-1,  3.845707852011820e-3, 4.112633446959460e-2,
         6.594735004847470e-1,
        ])

        E = np.array([
        -4.364293924074990e-2, -9.449231893318330e-02, -2.046776012570780e-2,
        -8.241166396760410e-2, -7.746685789572230e-02,  3.593817905677970e-2,
        ])

        F = np.array([-1.02956203215198])

        Eeq = A[0] * np.tanh((x + A[1]) / A[2]) \
            + A[3] * np.tanh((x + A[4]) / A[5]) \
            + A[6] * np.tanh((x + A[7]) / A[8]) \
            + A[9] * np.tanh((x + A[10]) / A[11]) \
            + A[12] * np.tanh((x + A[13]) / A[14]) \
            + A[15] * np.tanh((x + A[16]) / A[17]) \
            + A[18] \
            + B[0] * np.tanh((x + B[1]) / B[2]) \
            + B[3] * np.tanh((x + B[4]) / B[5]) \
            + (C[0] * x**8 + C[1] * x**7 + C[2] * x**6 + C[3] * x**5
               + C[4] * x**4 + C[5] * x**3 + C[6] * x**2 + C[7] * x + C[8]
               + D[0] * np.tanh((x + D[1]) / D[2])
               + D[3] * np.tanh((x + D[4]) / D[5])
               + D[6] * np.tanh((x + D[7]) / D[8])
               + D[9] * np.tanh((x + D[10]) / D[11])
               + D[12] * np.tanh((x + D[13]) / D[14])
               + D[15] * np.tanh((x + D[16]) / D[17])
               + D[18]
               + E[0] * np.tanh((x + E[1]) / E[2])
               + E[3] * np.tanh((x + E[4]) / E[5])) \
            / (1.0 + np.exp(-1.0e2 * (x + F[0])))

        Eeq = np.real(Eeq)

        Eeq[x <= 0] = 10
        Eeq[x > 1] = -10

        return Eeq.squeeze()

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

        M_hyst = 0.03
        if isinstance(x, np.ndarray):
            M_hyst *= np.ones_like(x)

        return M_hyst


class GraphiteSlow(GraphiteFast):

    def __init__(self, alpha_a: float, alpha_c: float, Li_max: float) -> None:
        """
        Computationally fast Graphite kinetic and transport properties.

        Differs from ``GraphiteSlow`` because the equilibrium potential is
        not piecewise here, making it less accurate, but faster to evaluate.

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

        super().__init__(alpha_a, alpha_c, Li_max)

        csvfile = os.path.dirname(__file__) + '/data/graphite_ocv.csv'
        df = pd.read_csv(csvfile).sort_values(by='x')

        self.x_min = df['x'].min()
        self.x_max = df['x'].max()
        self._Eeq_spline = CubicSpline(df['x'], df['V'])

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

        if isinstance(x, float):
            if x < self.x_min or x > self.x_max:
                raise ValueError(f'x is out of bounds [{self.x_min},'
                                 f' {self.x_max}]')

        if not isinstance(x, float):
            if np.any(x < self.x_min) or np.any(x > self.x_max):
                raise ValueError(f'x is out of bounds [{self.x_min},'
                                 f' {self.x_max}]')

        return self._Eeq_spline(x)


class GraphiteSlowExtrap(GraphiteFast):

    def __init__(self, alpha_a: float, alpha_c: float, Li_max: float) -> None:
        """
        Computationally fast Graphite kinetic and transport properties.

        Differs from ``GraphiteSlow`` because the piecewise equilibrium
        potential is extrapolated to be valid on the full [0, 1] range.

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

        super().__init__(alpha_a, alpha_c, Li_max)

        csvfile = os.path.dirname(__file__) + '/data/graphite_ocv_extrap.csv'
        df = pd.read_csv(csvfile).sort_values(by='x')

        self.x_min = df['x'].min()
        self.x_max = df['x'].max()
        self._Eeq_spline = CubicSpline(df['x'], df['V'])

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
        return self._Eeq_spline(x)
