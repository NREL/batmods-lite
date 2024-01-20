from numpy import ndarray as _ndarray


class NMC532Fast(object):

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

    def get_Ds(self, x: float | _ndarray, T: float) -> float | _ndarray:
        """
        Calculate the lithium diffusivity in the solid phase given the local
        intercalation fraction ``x`` and temperature ``T``.

        Parameters
        ----------
        x : float | 1D array
            Lithium intercalation fraction [-].

        T : float
            Battery temperature [K].

        Returns
        -------
        Ds : float | 1D array
            Lithium diffusivity in the solid phase [m^2/s].
        """

        import numpy as np

        from .. import Constants

        c = Constants()

        A = np.array([
        -2.509010843479270e+2,  2.391026725259970e+3, -4.868420267611360e+3,
        -8.331104102921070e+1,  1.057636028329000e+4, -1.268324548348120e+4,
         5.016272167775530e+3,  9.824896659649480e+2, -1.502439339070900e+3,
         4.723709304247700e+2, -6.526092046397090e+1])

        Ds = np.exp(-30e6 / c.R * (1 / T - 1 / 303.15)) \
           * 2.25 * 10.0**(  A[0] * x**10 + A[1] * x**9 + A[2] * x**8
                           + A[3] * x**7 + A[4] * x**6 + A[5] * x**5
                           + A[6] * x**4 + A[7] * x**3 + A[8] * x**2
                           + A[9] * x + A[10]  )

        return Ds

    def get_i0(self, x: float | _ndarray, C_Li: float | _ndarray,
               T: float) -> float | _ndarray:
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

        Returns
        -------
        i0 : float | 1D array
            Exchange current density [A/m^2].
        """

        import numpy as np

        from .. import Constants

        c = Constants()

        A = np.array([
         1.650452829641290e+1, -7.523567141488800e+1,  1.240524690073040e+2,
        -9.416571081287610e+1,  3.249768821737960e+1, -3.585290065824760e+0])

        i0 = 9 * (C_Li / 1.2)**self.alpha_a \
           * np.exp(-30e6 / c.R * (1 / T - 1 / 303.15)) \
           * (  A[0] * x**5 + A[1] * x**4 + A[2] * x**3 + A[3] * x**2
              + A[4] * x**1 + A[5]  )

        return i0

    def get_Eeq(self, x: float | _ndarray, T: float) -> float | _ndarray:
        """
        Calculate the equilibrium potential given the intercalation
        fraction ``x`` at the particle surface and temperature ``T``.

        Parameters
        ----------
        x : float | 1D array
            Lithium intercalation fraction at ``r = R_s`` [-].

        T : float
            Battery temperature [K].

        Returns
        -------
        Eeq : float | 1D array
            Equilibrium potential [V].
        """

        import numpy as np

        A = np.array([
         5.314735633000300e+0, -3.640117692001490e+3,  1.317657544484270e+4,
        -1.455742062291360e+4, -1.571094264365090e+3,  1.265630978512400e+4,
        -2.057808873526350e+3, -1.074374333186190e+4,  8.698112755348720e+3,
        -8.297904604107030e+2, -2.073765547574810e+3,  1.190223421193310e+3,
        -2.724851668445780e+2,  2.723409218042130e+1, -4.158276603609060e+0,
        -5.573191762723310e-4,  6.560240842659690e+0,  4.148209275061330e+1])

        Eeq = A[0] * x**0 + A[1] * x**14 + A[2] * x**13 + A[3] * x**12 \
            + A[4] * x**11 + A[5] * x**10 + A[6] * x**9 + A[7] * x**8  \
            + A[8] * x**7 + A[9] * x**6 + A[10] * x**5 + A[11] * x**4 \
            + A[12] * x**3 + A[13] * x**2 + A[14] * x**1 \
            + A[15] * np.exp( A[16] * x**A[17] )

        return Eeq
