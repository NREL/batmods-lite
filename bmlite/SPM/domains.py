"""
Domains Module
--------------
Contains classes to construct the battery for P2D simulations. Each class reads
in keyword arguments that define parameters relevant to its specific domain.
For example, the area and temperature are ``Battery`` level parameters because
they are the same everywhere, but the discretizations ``Nx`` and ``Nr`` may be
different for the anode, separator, and cathode domains.
"""

from numpy import ndarray as _ndarray


class Battery(object):

    def __init__(self, **kwargs) -> None:
        """
        A class for battery-level attributes.

        Parameters
        ----------
        **kwargs : dict, required
            Keyword arguments to set the battery-level attributes. The required
            keys and descriptions are given below:

            ====== =====================================================
            Key    Description [units] (type)
            ====== =====================================================
            cap    nominal battery capacity [A*h] (*float*)
            temp   temperature [K] (*float*)
            area   area normal to the current collectors [m^2] (*float*)
            ====== =====================================================
        """

        self.cap = kwargs.get('cap')
        self.temp = kwargs.get('temp')
        self.area = kwargs.get('area')

        self.update()

    def update(self):
        """
        Updates any secondary/dependent parameters. At present, this method
        does not do anything for the ``Battery`` class.

        Returns
        -------
        None.
        """
        pass


class Electrolyte(object):

    def __init__(self, **kwargs) -> None:
        """
        A class for the electrolyte attributes and methods.

        Parameters
        ----------
        **kwargs : dict, required
            Keyword arguments to set the electrolyte attributes. The required
            keys and descriptions are given below:

            ======== ==============================================
            Key      Description [units] (type)
            ======== ==============================================
            li_0     initial Li+ concentration [kmol/m^3] (*float*)
            ======== ==============================================
        """

        self.Li_0 = kwargs.get('Li_0', 1.2)

    def update(self) -> None:
        """
        Updates any secondary/dependent parameters. At present, this method
        does not do anything for the ``Electrolyte`` class.

        Returns
        -------
        None.
        """
        pass

    def make_mesh(self, pshift: int = 0) -> None:
        # [[ptr_an], [ptr_ca], phi_el]

        self.ptr = {}
        self.ptr['phi_el'] = 0 + pshift

        self.ptr['shift'] = 1

    def sv_0(self) -> _ndarray:
        import numpy as np
        return np.array([self.phi_0])

    def algidx(self) -> _ndarray:
        import numpy as np
        return np.hstack([self.ptr['phi_el']])

    def to_dict(self, sol) -> dict:

        el_sol = {}
        el_sol['phie'] = sol.y[:, self.ptr['phi_el']]

        return el_sol


class Electrode(object):

    def __init__(self, **kwargs):
        """
        A class for the electrode-specific attributes and methods.

        This class is used to build both the anode and cathode for P2D
        simulations.

        Parameters
        ----------
        **kwargs : dict, required
            Keyword arguments to set the electrode attributes. The required
            keys and descriptions are given below:

            ========= =========================================================
            Key       Description [units] (type)
            ========= =========================================================
            Nr        number of ``r`` discretizations [-] (*int*)
            thick     electrode thickness [m] (*float*)
            R_s       represenatative particle radius [m] (*float*)
            eps_el    electrolyte volume fraction [-] (*float*)
            eps_CBD   carbon binder domain volume fraction [-] (*float*)
            eps_void  Void volume fraction [-] (*float*)
            alpha_a   Butler-Volmer anodic symmetry factor [-] (*float*)
            alpha_c   Butler-Volmer cathodic symmetry factor [-] (*float*)
            Li_max    max solid-phase lithium concentraion [kmol/m^3] (*float*)
            x_0       initial solid-phase intercalation fraction [-] (*float*)
            i0_deg    ``i0`` degradation factor [-] (*float*)
            Ds_deg    ``Ds`` degradation factor [-] (*float*)
            material  class name from ``bmlite.materials`` [-] (*str*)
            ========= =========================================================
        """

        self.Nr = kwargs.get('Nr')
        self.thick = kwargs.get('thick')
        self.R_s = kwargs.get('R_s')
        self.eps_el = kwargs.get('eps_el')
        self.eps_CBD = kwargs.get('eps_CBD')
        self.eps_void = kwargs.get('eps_void')
        self.alpha_a = kwargs.get('alpha_a')
        self.alpha_c = kwargs.get('alpha_c')
        self.Li_max = kwargs.get('Li_max')
        self.x_0 = kwargs.get('x_0')
        self.i0_deg = kwargs.get('i0_deg')
        self.Ds_deg = kwargs.get('Ds_deg')
        self.material = kwargs.get('material')

        self.update()

    def update(self) -> None:
        """
        Updates any secondary/dependent parameters. For the ``Electrode``
        class, this initializes the material class, and sets the following:

        * Solid-phase volume fraction [-]:
            ``eps_s = 1 - eps_el``
        * Activate material volume fraction [-]:
            ``eps_AM = 1 - eps_el - eps_void - eps_CBD``
        * Specific particle surface area [m^2/m^3]:
            ``A_s = eps_AM / R_s``

        Returns
        -------
        None.
        """

        from .. import materials

        self.eps_s = 1 - self.eps_el
        self.eps_AM = 1 - self.eps_el - self.eps_void - self.eps_CBD
        self.A_s = 3 * self.eps_AM / self.R_s

        Material = getattr(materials, self.material)
        self._material = Material(self.alpha_a, self.alpha_c, self.Li_max)

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

        return self.Ds_deg * self._material.get_Ds(x, T)

    def get_i0(self, x: float, C_Li: float, T: float) -> float:
        """
        Calculate the exchange current density given the surface intercalation
        fraction ``x`` at the particle surface, the local lithium ion
        concentration ``C_Li``, and temperature ``T``.

        Parameters
        ----------
        x : float
            Lithium intercalation fraction at ``r = R_s`` [-].

        C_Li : float
            Lithium ion concentration in the local electrolyte [kmol/m^3].

        T : float
            Battery temperature [K].

        Returns
        -------
        i0 : float
            Exchange current density [A/m^2].
        """

        return self.i0_deg * self._material.get_i0(x, C_Li, T)

    def get_Eeq(self, x: float, T: float) -> float:
        """
        Calculate the equilibrium potential given the surface intercalation
        fraction ``x`` at the particle surface and temperature ``T``.

        Parameters
        ----------
        x : float
            Lithium intercalation fraction at ``r = R_s`` [-].

        T : float
            Battery temperature [K].

        Returns
        -------
        Eeq : float
            Equilibrium potential [V].
        """

        return self._material.get_Eeq(x, T)

    def make_mesh(self, pshift: int = 0):
        """
        Determines/sets the ``r`` locations for all of the "minus" interfaces
        ``rm``, "plus" interfaces ``rp``, and control volume centers ``r``
        based on the representative particle radius and ``Nr`` discretization.
        At present, only a uniform mesh is supported.

        Parameters
        ----------
        pshift : int, optional

        Returns
        -------
        None.

        See also
        --------
        batmods.mesh.r_ptr, batmods.mesh.uniform_mesh
        """

        from ..mesh import r_ptr, uniform_mesh

        # Mesh locations
        self.rm, self.rp, self.r = uniform_mesh(self.R_s, self.Nr)

        # Pointers
        # [[ptr_an], [ptr_ca], phi_el]
        # ptr_an and ptr_ca -> [[Li_ed(0->R_s)], phi_ed]

        self.ptr = {}
        self.ptr['Li_ed'] = 0 + pshift
        self.ptr['r_off'] = 1

        self.ptr['phi_ed'] = self.ptr['Li_ed'] + self.Nr

        self.ptr['shift'] = self.Nr + 1

        r_ptr(self, ['Li_ed'])

    def sv_0(self):
        import numpy as np
        return np.hstack([self.x_0 * np.ones(self.Nr), self.phi_0])

    def algidx(self):
        import numpy as np
        return np.array([self.ptr['phi_ed']])

    def to_dict(self, sol: object) -> dict:

        ed_sol = {}
        ed_sol['phis'] = sol.y[:, self.ptr['phi_ed']]
        ed_sol['xs'] = sol.y[:, self.r_ptr['Li_ed']]
        ed_sol['cs'] = ed_sol['xs'] * self.Li_max

        return ed_sol
