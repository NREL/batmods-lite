"""
Domains Module
--------------
Contains classes to construct the battery for SPM simulations. Each class reads
in keyword arguments that define parameters relevant to its specific domain.
For example, the area and temperature are ``Battery`` level parameters because
they are the same everywhere, but the discretization ``Nr`` may be different
for the anode and cathode domains.

"""

import numpy as np


class Battery:

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
            area   area normal to the current collectors [m2] (*float*)
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


class Electrolyte:

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
            li_0     initial Li+ concentration [kmol/m3] (*float*)
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
        self.ptr['phie'] = 0 + pshift

        self.ptr['shift'] = 1

    def sv0(self) -> np.ndarray:
        return np.array([self.phi_0])

    def algidx(self) -> np.ndarray:
        return np.hstack([self.ptr['phie']])

    def to_dict(self, soln) -> dict:

        el_soln = {}
        el_soln['phie'] = soln.y[:, self.ptr['phie']]

        return el_soln


class Electrode:

    def __init__(self, name: str, **kwargs):
        """
        A class for the electrode-specific attributes and methods.

        This class is used to build both the anode and cathode for P2D
        simulations.

        Parameters
        ----------
        name : str
            Name of the electrode, must be either 'anode' or 'cathode'.
        **kwargs : dict, required
            Keyword arguments to set the electrode attributes. The required
            keys and descriptions are given below:

            ========== ======================================================
            Key        Description [units] (type)
            ========== ======================================================
            Nr         number of ``r`` discretizations [-] (*int*)
            thick      electrode thickness [m] (*float*)
            R_s        represenatative particle radius [m] (*float*)
            eps_el     electrolyte volume fraction [-] (*float*)
            eps_CBD    carbon binder domain volume fraction [-] (*float*)
            eps_void   Void volume fraction [-] (*float*)
            alpha_a    Butler-Volmer anodic symmetry factor [-] (*float*)
            alpha_c    Butler-Volmer cathodic symmetry factor [-] (*float*)
            Li_max     max solid-phase Li concentration [kmol/m3] (*float*)
            x_0        initial solid Li intercalation fraction [-] (*float*)
            i0_deg     ``i0`` degradation factor [-] (*float*)
            Ds_deg     ``Ds`` degradation factor [-] (*float*)
            material   class name from ``bmlite.materials`` [-] (*str*)
            submodels  ``submodels`` classes to include (*dict[dict]*)
            ========== ======================================================

        """

        from . import submodels

        if name not in ['anode', 'cathode']:
            raise ValueError("'name' must be either 'anode' or 'cathode'.")

        self._name = name

        self.Nr = kwargs.get('Nr')
        self.thick = kwargs.get('thick')
        self.R_s = kwargs.get('R_s')
        self.eps_s = kwargs.get('eps_s')
        self.eps_el = kwargs.get('eps_el')
        self.eps_CBD = kwargs.get('eps_CBD')
        self.alpha_a = kwargs.get('alpha_a')
        self.alpha_c = kwargs.get('alpha_c')
        self.Li_max = kwargs.get('Li_max')
        self.x_0 = kwargs.get('x_0')
        self.i0_deg = kwargs.get('i0_deg')
        self.Ds_deg = kwargs.get('Ds_deg')
        self.material = kwargs.get('material')

        self.update()

        self._submodels = {}
        all_submodels = kwargs.get('submodels', {})
        if 'Hysteresis' in all_submodels:
            Hysteresis, opt = submodels.Hysteresis, all_submodels['Hysteresis']
            self._submodels['Hysteresis'] = Hysteresis(self, **opt)

    def update(self) -> None:
        """
        Updates any secondary/dependent parameters. For the ``Electrode``
        class, this initializes the material class, and sets the following:

        * Void-phase volume fraction [-]:
            ``eps_void = 1 - eps_s - eps_el``
        * Activate material volume fraction [-]:
            ``eps_AM = eps_s - eps_CBD``
        * Specific particle surface area [m2/m3]:
            ``A_s = 3 * eps_AM / R_s``

        Returns
        -------
        None.

        """

        from .. import materials

        self.eps_void = 1. - self.eps_s - self.eps_el
        self.eps_AM = self.eps_s - self.eps_CBD
        self.A_s = 3. * self.eps_AM / self.R_s

        if self.eps_void < 0.:
            raise ValueError('eps_s + eps_el > 1.0')

        Material = getattr(materials, self.material)
        self._material = Material(self.alpha_a, self.alpha_c, self.Li_max)

    def get_Ds(self, x: float | np.ndarray, T: float,
               fluxdir: float) -> float | np.ndarray:
        """
        Calculate the lithium diffusivity in the solid phase given the local
        intercalation fraction ``x`` and temperature ``T``.

        Parameters
        ----------
        x : float | 1D array
            Lithium intercalation fraction [-].
        T : float
            Battery temperature [K].
        fluxdir : float
            Lithiation direction: +1 for lithiation, -1 for delithiation, 0 for
            rest. Used for direction-dependent parameters. Ensure the zero case
            is handled explicitly or via a default (lithiating or delithiating).

        Returns
        -------
        Ds : float | 1D array
            Lithium diffusivity in the solid phase [m2/s].

        """
        return self.Ds_deg * self._material.get_Ds(x, T, fluxdir)

    def get_i0(self, x: float, C_Li: float, T: float,
               fluxdir: float) -> float:
        """
        Calculate the exchange current density given the surface intercalation
        fraction ``x`` at the particle surface, the local lithium ion
        concentration ``C_Li``, and temperature ``T``.

        Parameters
        ----------
        x : float
            Lithium intercalation fraction at ``r = R_s`` [-].
        C_Li : float
            Lithium ion concentration in the local electrolyte [kmol/m3].
        T : float
            Battery temperature [K].
        fluxdir : float
            Lithiation direction: +1 for lithiation, -1 for delithiation, 0 for
            rest. Used for direction-dependent parameters. Ensure the zero case
            is handled explicitly or via a default (lithiating or delithiating).

        Returns
        -------
        i0 : float
            Exchange current density [A/m2].

        """
        return self.i0_deg * self._material.get_i0(x, C_Li, T, fluxdir)

    def get_Eeq(self, x: float) -> float:
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
        return self._material.get_Eeq(x)

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
        return self._material.get_Mhyst(x)

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

        self._wtm = 0.5*(self.rp[:-1] - self.rm[:-1]) / np.diff(self.r)

        self._wtp = 0.5*(self.rp[1:] - self.rm[1:]) / np.diff(self.r)

        # Pointers
        # [[ptr_an], phi_el, [ptr_ca]]
        # ptr_an -> [[Li_ed(0->R_s)], phi_ed, ...]
        # ptr_ca -> [..., phi_ed, [Li_ed(R_s->0)]]

        self.ptr = {}
        if self._name == 'anode':
            self.ptr['xs'] = 0 + pshift
            self.ptr['phis'] = self.ptr['xs'] + self.Nr

        elif self._name == 'cathode':
            self.ptr['phis'] = 0 + pshift
            self.ptr['xs'] = self.ptr['phis'] + 1

        for model in self._submodels.values():
            model.make_mesh(pshift)

        submodel_count = len(self._submodels)

        self.ptr['r_off'] = 1
        self.ptr['start'] = pshift
        self.ptr['size'] = self.Nr + 1 + submodel_count
        self.ptr['shift'] = self.ptr['size']

        r_ptr(self, ['xs'])

    def sv0(self):

        start = self.ptr['start']
        size = self.ptr['size']

        sv0 = np.zeros(size)
        sv0[self.r_ptr['xs'] - start] = self.x_0*np.ones(self.Nr)
        sv0[self.ptr['phis'] - start] = self.phi_0

        for model in self._submodels.values():
            model.sv0(sv0)

        return sv0

    def algidx(self):

        algidx = np.array([self.ptr['phis']], dtype=int)
        for model in self._submodels.values():
            model.algidx(algidx)

        return np.sort(algidx)

    def to_dict(self, soln: object) -> dict:

        phis = soln.y[:, self.ptr['phis']]

        xs = soln.y[:, self.r_ptr['xs']]
        if self._name == 'cathode':
            xs = np.flip(xs, axis=1)

        ed_soln = {
            'r': self.r,
            'phis': phis,
            'xs': xs,
            'cs': xs*self.Li_max,
        }

        for model in self._submodels.values():
            outputs = model.to_dict(soln)
            ed_soln.update(outputs)

        return ed_soln

    def _boundary_voltage(self, soln) -> np.ndarray:
        """
        Calculate and return the boundary voltage at all solution times.

        Parameters
        ----------
        soln : Solution
            A SPM solution instance, step or cycle.

        Returns
        -------
        voltage_V : np.ndarray
            Boundary voltage, in volts.

        """

        V_ptr = self.ptr['phis']

        return soln.y[:, V_ptr]

    def _boundary_current(self, soln) -> np.ndarray:
        """
        Calculate and return the boundary current at all solution times.

        Parameters
        ----------
        soln : Solution
            A SPM solution instance, step or cycle.

        Returns
        -------
        current_A : np.ndarray
            Boundary current, in amps.

        """

        sim = soln._sim

        bat, el = sim.bat, sim.el

        c, T = sim.c, bat.temp

        # calculate the boundary current using sum of Fardaic reactions
        if self._name == 'anode':
            sign = +1.
        elif self._name == 'cathode':
            sign = -1.
        else:
            raise ValueError("Electrode name not in {'anode', 'cathode'}.")

        ed_soln = self.to_dict(soln)
        el_soln = el.to_dict(soln)

        phis = ed_soln['phis']
        xs_R = ed_soln['xs'][:, -1]
        phie = el_soln['phie']

        if 'Hysteresis' in self._submodels:
            hyst = ed_soln['hyst']
            Hyst = self.get_Mhyst(xs_R)*hyst
        else:
            Hyst = 0.

        eta = phis - phie - (self.get_Eeq(xs_R) + Hyst)
        fluxdir = -np.sign(eta)

        i0 = self.get_i0(xs_R, el.Li_0, T, fluxdir)
        sdot = i0 / c.F * (  np.exp( self.alpha_a*c.F*eta / c.R / T)
                           - np.exp(-self.alpha_c*c.F*eta / c.R / T)  )

        i_ext = sign*sdot*self.A_s*self.thick*c.F
        current_A = i_ext*bat.area

        return current_A
