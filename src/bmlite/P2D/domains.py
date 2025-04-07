"""
Domains Module
--------------
Contains classes to construct the battery for P2D simulations. Each class reads
in keyword arguments that define parameters relevant to its specific domain.
For example, the area and temperature are ``Battery`` level parameters because
they are the same everywhere, but the discretizations ``Nx`` and ``Nr`` may be
different for the anode, separator, and cathode domains.

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

            ====== =================================================
            Key    Description [units] (type)
            ====== =================================================
            cap    nominal battery capacity [A*h] (*float*)
            temp   temperature [K] (*float*)
            area   area normal to current collectors [m^2] (*float*)
            ====== =================================================

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

            ========== ================================================
            Key        Description [units] (type)
            ========== ================================================
            li_0       initial Li+ concentration [kmol/m^3] (*float*)
            material   class name from ``bmlite.materials`` [-] (*str*)
            ========== ================================================

        """

        self.Li_0 = kwargs.get('Li_0')
        self.material = kwargs.get('material')

    def update(self) -> None:
        """
        Updates any secondary/dependent parameters. For the ``Electrolyte``
        class, this only initializes the material class.

        Returns
        -------
        None.

        """

        from .. import materials

        ElyteMaterial = getattr(materials, self.material)
        self._material = ElyteMaterial()

    def get_D(self, C_Li: float | np.ndarray, T: float) -> float | np.ndarray:
        """
        Calculate the lithium ion diffusivity in the electrolyte solution at
        concentration ``C_Li`` and temperature ``T``.

        Parameters
        ----------
        C_Li : float | 1D array
            Lithium ion concentration in the electrolyte [kmol/m^3].
        T : float
            Battery temperature [K].

        Returns
        -------
        D : float | 1D array
            Lithium ion diffusivity in the electrolyte [m^2/s].

        """

        return self._material.get_D(C_Li, T)

    def get_t0(self, C_Li: float | np.ndarray, T: float) -> float | np.ndarray:
        """
        Calculate the lithium ion transference number at concentration ``C_Li``
        and temperature ``T``.

        Parameters
        ----------
        C_Li : float | 1D array
            Lithium ion concentration in the electrolyte [kmol/m^3].
        T : float
            Battery temperature [K].

        Returns
        -------
        t0 : float | 1D array
            Lithium ion transference number [-].

        """

        return self._material.get_t0(C_Li, T)

    def get_kappa(self, C_Li: float | np.ndarray,
                  T: float) -> float | np.ndarray:
        """
        Calculate the electrolyte conductivity at concentration ``C_Li`` and
        temperature ``T``.

        Parameters
        ----------
        C_Li : float | 1D array
            Lithium ion concentration in the electrolyte [kmol/m^3].
        T : float
            Battery temperature [K].

        Returns
        -------
        kappa : float | 1D array
            Electrolyte conductivity [S/m].

        """

        return self._material.get_kappa(C_Li, T)

    def get_gamma(self, C_Li: float | np.ndarray,
                  T: float) -> float | np.ndarray:
        """
        Calculate the electrolyte thermodynamic factor at concentration
        ``C_Li`` and temperature ``T``.

        Parameters
        ----------
        C_Li : float | 1D array
            Lithium ion concentration in the electrolyte [kmol/m^3].
        T : float
            Battery temperature [K].

        Returns
        -------
        gamma : float | 1D array
            Thermodynamic factor [-].

        """

        return self._material.get_gamma(C_Li, T)


class Electrode:

    def __init__(self, name: str, **kwargs) -> None:
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

            ========= ========================================================
            Key       Description [units] (type)
            ========= ========================================================
            Nx        number of ``x`` discretizations [-] (*int*)
            Nr        number of ``r`` discretizations [-] (*int*)
            thick     electrode thickness [m] (*float*)
            R_s       represenatative particle radius [m] (*float*)
            eps_el    electrolyte volume fraction [-] (*float*)
            eps_CBD   carbon binder domain volume fraction [-] (*float*)
            eps_void  void volume fraction [-] (*float*)
            p_sol     solid Bruggeman factor, ``eps_s**p_sol`` [-] (*float*)
            p_liq     liquid Bruggeman factor, ``eps_el**p_liq`` [-] (*float*)
            alpha_a   Butler-Volmer anodic symmetry factor [-] (*float*)
            alpha_c   Butler-Volmer cathodic symmetry factor [-] (*float*)
            Li_max    max solid-phase Li concentration [kmol/m^3] (*float*)
            x_0       initial solid-phase intercalation [-] (*float*)
            i0_deg    ``i0`` degradation factor [-] (*float*)
            Ds_deg    ``Ds`` degradation factor [-] (*float*)
            material  class name from ``bmlite.materials`` [-] (*str*)
            ========= ========================================================

        """

        from . import options

        if name not in ['anode', 'cathode']:
            raise ValueError("'name' must be either 'anode' or 'cathode'.")

        self.name = name

        self.Nx = kwargs.get('Nx')
        self.Nr = kwargs.get('Nr')
        self.thick = kwargs.get('thick')
        self.R_s = kwargs.get('R_s')
        self.eps_s = kwargs.get('eps_s')
        self.eps_el = kwargs.get('eps_el')
        self.eps_CBD = kwargs.get('eps_CBD')
        self.p_sol = kwargs.get('p_sol')
        self.p_liq = kwargs.get('p_liq')
        self.alpha_a = kwargs.get('alpha_a')
        self.alpha_c = kwargs.get('alpha_c')
        self.Li_max = kwargs.get('Li_max')
        self.x_0 = kwargs.get('x_0')
        self.i0_deg = kwargs.get('i0_deg')
        self.Ds_deg = kwargs.get('Ds_deg')
        self.material = kwargs.get('material')

        self._options = {}

        all_options = kwargs.get('options', {})
        if 'Hysteresis' in all_options:
            Hysteresis, opt = options.Hysteresis, all_options['Hysteresis']
            self._options['Hysteresis'] = Hysteresis(self, **opt)

        self.update()

    def update(self) -> None:
        """
        Updates any secondary/dependent parameters. For the ``Electrode``
        class, this initializes the material class, and sets the following:

        * Void-phase volume fraction [-]:
            ``eps_void = 1 - eps_s - eps_el``
        * Activate material volume fraction [-]:
            ``eps_AM = eps_s - eps_CBD``
        * Solid-phase conductivity [S/m]:
            ``sigma_s = 10 * eps_s``
        * Specific particle surface area [m^2/m^3]:
            ``A_s = 3 * eps_AM / R_s``

        Returns
        -------
        None.

        """

        from .. import materials

        self.eps_void = 1. - self.eps_s - self.eps_el
        self.eps_AM = self.eps_s - self.eps_CBD
        self.sigma_s = 10. * self.eps_s
        self.A_s = 3. * self.eps_AM / self.R_s

        if self.eps_void < 0.:
            raise ValueError('eps_s + eps_el > 1.0')

        ActiveMaterial = getattr(materials, self.material)
        self._material = ActiveMaterial(self.alpha_a, self.alpha_c,
                                        self.Li_max)

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

        return self.Ds_deg * self._material.get_Ds(x, T, fluxdir)

    def get_i0(self, x: float | np.ndarray, C_Li: float | np.ndarray,
               T: float, fluxdir: float | np.ndarray) -> float | np.ndarray:
        """
        Calculate the exchange current density given the surface intercalation
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

        return self.i0_deg * self._material.get_i0(x, C_Li, T, fluxdir)

    def get_Eeq(self, x: float | np.ndarray) -> float | np.ndarray:
        """
        Calculate the equilibrium potential given the surface intercalation
        fraction ``x`` at the particle surface and temperature ``T``.

        Parameters
        ----------
        x : float | 1D array
            Lithium intercalation fraction at ``r = R_s`` [-].

        Returns
        -------
        Eeq : float | 1D array
            Equilibrium potential [V].

        """

        return self._material.get_Eeq(x)

    def make_mesh(self, xshift: float = 0., pshift: int = 0) -> None:
        """
        Determines/sets the ``x`` locations for all of the "minus" interfaces
        ``xm``, "plus" interfaces ``xp``, and control volume centers ``x``
        based on the electrode thickness and ``Nx`` discretization. At present,
        only a uniform mesh is supported.

        Also determines/sets the ``r`` locations for all of the "minus"
        interfaces ``rm``, "plus" interfaces ``rp``, and control volume centers
        ``r`` based on the representative particle radius and ``Nr``
        discretization. At present, only a uniform mesh is supported.

        Parameters
        ----------
        xshift : float, optional
        pshift : int, optional

        Returns
        -------
        None.

        See also
        --------
        batmods.mesh.x_ptr, batmods.mesh.xr_ptr, batmods.mesh.uniform_mesh

        """

        from ..mesh import x_ptr, xr_ptr, uniform_mesh

        # Mesh locations
        self.xm, self.xp, self.x = uniform_mesh(self.thick, self.Nx, xshift)
        self.rm, self.rp, self.r = uniform_mesh(self.R_s, self.Nr)

        # Pointers
        # [[ptr_an], [ptr_sep], [ptr_ca]]
        # ptr_an and ptr_ca -> [[Li_ed(0->R_s)], phi_ed, Li_el, phi_el, ...]

        self.ptr = {}
        self.ptr['Li_ed'] = 0 + pshift
        self.ptr['r_off'] = 1

        self.ptr['phi_ed'] = self.ptr['Li_ed'] + self.Nr
        self.ptr['Li_el'] = self.ptr['phi_ed'] + 1
        self.ptr['phi_el'] = self.ptr['Li_el'] + 1

        xvars = ['phi_ed', 'Li_el', 'phi_el']
        last_xvar = 'phi_el'

        # Options only support new x variables (like hysteresis... not xr)
        for opt in self._options.values():
            new_xvar = opt.make_mesh(last_xvar, pshift)
            xvars.append(new_xvar)
            last_xvar = new_xvar

        opt_count = len(self._options)

        self.ptr['x_off'] = self.Nr + 3 + opt_count

        self.ptr['start'] = pshift
        self.ptr['size'] = self.ptr['x_off']*self.Nx
        self.ptr['shift'] = self.ptr['size']

        x_ptr(self, xvars)
        xr_ptr(self, ['Li_ed'])

    def sv0(self, el: object) -> np.ndarray:

        start = self.ptr['start']
        size = self.ptr['size']

        sv0 = np.zeros(size)
        sv0[self.xr_ptr['Li_ed'].flatten() - start] = self.x_0
        sv0[self.x_ptr['phi_ed'] - start] = self.phi_0
        sv0[self.x_ptr['Li_el'] - start] = el.Li_0
        sv0[self.x_ptr['phi_el'] - start] = el.phi_0

        for opt in self._options.values():
            opt.sv0(sv0)

        return sv0

    def algidx(self) -> np.ndarray:
        algidx = np.hstack([self.x_ptr['phi_ed'], self.x_ptr['phi_el']])

        for opt in self._options.values():
            opt.algidx(algidx)

        return np.sort(algidx)

    def to_dict(self, sol: object) -> dict:

        xs = np.zeros([sol.t.size, self.Nx, self.Nr])
        for i in range(sol.t.size):
            xs_tmp = sol.y[i, self.xr_ptr['Li_ed'].flatten()]
            xs[i, :, :] = xs_tmp.reshape([self.Nx, self.Nr])

        ed_sol = {
            'x': self.x,
            'r': self.r,
            'xs': xs,
            'cs': xs*self.Li_max,
            'phis': sol.y[:, self.x_ptr['phi_ed']],
            'ce': sol.y[:, self.x_ptr['Li_el']],
            'phie': sol.y[:, self.x_ptr['phi_el']],
        }

        for opt in self._options.values():
            outputs = opt.to_dict(sol)
            ed_sol.update(outputs)

        return ed_sol

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

        bat, an, ca = sim.bat, sim.an, sim.ca

        c, T = sim.c, bat.temp

        if self.name == 'anode':
            ed, ed_str = an, 'an'
        elif self.name == 'cathode':
            ed, ed_str = ca, 'ca'

        phis = soln.vars[ed_str]['phis']
        xs = soln.vars[ed_str]['xs'][:, :, -1]
        phie = soln.vars[ed_str]['phie']
        ce = soln.vars[ed_str]['ce']

        if 'Hysteresis' in ed._options:
            hyst = soln.vars[ed_str]['hyst']
            Hyst = ed.M_hyst*hyst
        else:
            Hyst = 0.

        eta = phis - phie - (ed.get_Eeq(xs) + Hyst)
        fluxdir = -np.sign(eta)

        i0 = ed.get_i0(xs, ce, T, fluxdir)
        sdot = i0 / c.F * (  np.exp( ed.alpha_a*c.F*eta / c.R / T)
                           - np.exp(-ed.alpha_c*c.F*eta / c.R / T)  )

        if self.name == 'anode':
            i_ext = sdot[:, 0]*an.A_s*c.F*(an.xp[0] - an.xm[0]) \
                    - an.sigma_s*an.eps_s**an.p_sol \
                        * (phis[:, 1] - phis[:, 0]) / (an.x[1] - an.x[0])

        elif self.name == 'cathode':
            i_ext = -sdot[:, -1]*ca.A_s*c.F*(ca.xp[-1] - ca.xm[-1]) \
                  - ca.sigma_s*ca.eps_s**ca.p_sol \
                      * (phis[:, -1] - phis[:, -2]) / (ca.x[-1] - ca.x[-2])

        current_A = i_ext*bat.area

        return current_A


class Separator:

    def __init__(self, **kwargs) -> None:
        """
        A class for the separator attributes and methods.

        This class is used to build both the separator for P2D simulations.

        Parameters
        ----------
        **kwargs : dict, required
            Keyword arguments to set the separator attributes. The required
            keys and descriptions are given below:

            ======== ========================================================
            Key      Description [units] (type)
            ======== ========================================================
            Nx       number of ``x`` discretizations [-] (*int*)
            thick    electrode thickness [m] (*float*)
            eps_el   electrolyte volume fraction [-] (*float*)
            p_liq    liquid Bruggeman factor, ``eps_el**p_liq`` [-] (*float*)
            ======== ========================================================

        """

        self.Nx = kwargs.get('Nx')
        self.thick = kwargs.get('thick')
        self.eps_el = kwargs.get('eps_el')
        self.p_liq = kwargs.get('p_liq')

        self.update()

    def update(self) -> None:
        """
        Updates any secondary/dependent parameters. For the ``Separator``
        class, this sets the following:

        * Solid-phase volume fraction [-]:
            ``eps_s = 1 - eps_el``

        Returns
        -------
        None.

        """

        self.eps_s = 1 - self.eps_el

    def make_mesh(self, xshift: float = 0., pshift: int = 0) -> None:
        """
        Determines/sets the ``x`` locations for all of the "minus" interfaces
        ``xm``, "plus" interfaces ``xp``, and control volume centers ``x``
        based on the electrode thickness and ``Nx`` discretization. At present,
        only a uniform mesh is supported.

        Parameters
        ----------
        xshift : float, optional
        pshift : int, optional

        Returns
        -------
        None.

        See also
        --------
        batmods.mesh.x_ptr, batmods.mesh.xr_ptr, batmods.mesh.uniform_mesh

        """

        from ..mesh import x_ptr, uniform_mesh

        # Mesh locations
        self.xm, self.xp, self.x = uniform_mesh(self.thick, self.Nx, xshift)

        # Pointers
        # [[ptr_an], Li_el, phi_el, ..., [ptr_ca]]

        self.ptr = {}
        self.ptr['Li_el'] = 0 + pshift
        self.ptr['phi_el'] = self.ptr['Li_el'] + 1
        self.ptr['x_off'] = 2

        self.ptr['shift'] = self.Nx * self.ptr['x_off']

        x_ptr(self, ['Li_el', 'phi_el'])

    def sv0(self, el: object) -> np.ndarray:
        import numpy as np
        return np.tile([el.Li_0, el.phi_0], self.Nx)

    def algidx(self) -> np.ndarray:
        return self.x_ptr['phi_el']

    def to_dict(self, sol: object) -> dict:

        sep_sol = {
            'x': self.x,
            'ce': sol.y[:, self.x_ptr['Li_el']],
            'phie': sol.y[:, self.x_ptr['phi_el']],
        }

        return sep_sol
