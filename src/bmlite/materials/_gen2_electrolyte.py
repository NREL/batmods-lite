from numpy import ndarray as _ndarray


class Gen2Electrolyte:

    def __init__(self) -> None:
        """
        Gen 2 Electrolyte transport and material properties.
        """
        pass

    def get_D(self, C_Li: float | _ndarray, T: float) -> float | _ndarray:
        """
        Calculate the lithium ion diffusivity in the electrolyte solution at
        concentration ``C_Li`` and temperature ``T``.

        Parameters
        ----------
        C_Li : float | 1D array
            Lithium ion concentration in the electrolyte [kmol/m3].

        T : float
            Battery temperature [K].

        Returns
        -------
        D : float | 1D array
            Lithium ion diffusivity in the electrolyte [m2/s].
        """

        import numpy as np

        A = np.array([[-0.568822600, 1607.003, -24.83763, 64.07366],
                      [-0.810872100, 475.2910, -24.83763, 64.07366],
                      [-0.005192312, 33.43827, -24.83763, 64.07366]])

        D = 0.0001 * 10**(
            (A[0, 0] - A[0, 1] / (T - (A[0, 2] + A[0, 3] * C_Li)))
            + (A[1, 0] + A[1, 1] / (T - (A[1, 2] + A[1, 3] * C_Li))) * C_Li
            + (A[2, 0] - A[2, 1] / (T - (A[2, 2] + A[2, 3] * C_Li))) * C_Li**2)

        return D

    def get_t0(self, C_Li: float | _ndarray, T: float) -> float | _ndarray:
        """
        Calculate the lithium ion transference number at concentration ``C_Li``
        and temperature ``T``.

        Parameters
        ----------
        C_Li : float | 1D array
            Lithium ion concentration in the electrolyte [kmol/m3].

        T : float
            Battery temperature [K].

        Returns
        -------
        t0 : float | 1D array
            Lithium ion transference number [-].
        """

        import numpy as np

        A = np.array([[-0.0000002876102,  0.0002077407, -0.03881203],
                      [ 0.0000011614630, -0.0008682500,  0.17772660],
                      [-0.0000006766258,  0.0006389189,  0.30917610]])

        t0 = np.polyval(A[0, :], T)*C_Li**2 + np.polyval(A[1, :], T)*C_Li \
           + np.polyval(A[2, :], T)

        return t0

    def get_kappa(self, C_Li: float | _ndarray, T: float) -> float | _ndarray:
        """
        Calculate the electrolyte conductivity at concentration ``C_Li`` and
        temperature ``T``.

        Parameters
        ----------
        C_Li : float | 1D array
            Lithium ion concentration in the electrolyte [kmol/m3].

        T : float
            Battery temperature [K].

        Returns
        -------
        kappa : float | 1D array
            Electrolyte conductivity [S/m].
        """

        import numpy as np

        A = np.array([
        [ 0.,           0.,           1.909446e-4, -8.038545e-2,  9.003410e+0],
        [-2.887587e-8,  3.483638e-5, -1.583677e-2,  3.195295e+0, -2.414638e+2],
        [ 1.653786e-8, -1.99876e-5,   9.071155e-3, -1.828064e+0,  1.380976e+2],
        [-2.791965e-9,  3.377143e-6, -1.532707e-3,  3.090003e-1, -2.335671e+1],
        ])

        kappa = np.polyval(A[0, :], T)*C_Li + np.polyval(A[1, :], T)*C_Li**2 \
              + np.polyval(A[2, :], T)*C_Li**3 + np.polyval(A[3, :], T)*C_Li**4

        return kappa

    def get_gamma(self, C_Li: float | _ndarray, T: float) -> float | _ndarray:
        """
        Calculate the electrolyte thermodynamic factor at concentration
        ``C_Li`` and temperature ``T``.

        Parameters
        ----------
        C_Li : float | 1D array
            Lithium ion concentration in the electrolyte [kmol/m3].

        T : float
            Battery temperature [K].

        Returns
        -------
        gamma : float | 1D array
            Thermodynamic factor [-].
        """

        import numpy as np

        gamma = 0.54000*np.exp(329./T)*C_Li**2 - 0.00225*np.exp(1360./T)*C_Li \
              + 0.34100*np.exp(261./T)

        return gamma
