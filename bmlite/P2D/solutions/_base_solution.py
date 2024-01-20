from numpy import ndarray as _ndarray


class BaseSolution(object):
    """
    Base methods for all P2D Solution classes.
    """

    __slots__ = ['_sim', '_exp', '_t', '_y', '_ydot', '_success', '_onroot',
                 '_message', '_solvetime']

    def __init__(self, sim: object, exp: dict) -> None:
        """
        When initialized, a copy of the P2D simulation object and experimental
        details for the solution instance are stored, all other instance
        attributes are set to ``None``. These are filled later with a "fill"
        method (e.g., ``ida_fill()`` or ``dict_fill()``).

        Parameters
        ----------
        sim : P2D Simulation object
            The P2D Simulation instance used to produce the solution.

        exp : dict
            Experiment dictionary. Specific key/value pairs are dependent on
            the experiment that was run.
        """

        self._sim = sim.copy()
        self._exp = exp.copy()

        self._t = None
        self._y = None
        self._ydot = None
        self._success = None
        self._onroot = None
        self._message = None
        self._solvetime = None

    @property
    def classname(self) -> str:
        """
        Class name.

        Returns
        -------
        classname : str
            Name of current class.
        """
        return 'BaseSolution'

    @property
    def t(self) -> _ndarray:
        """
        Saved solution times [s].

        Returns
        -------
        t : 1D array
            Solution times [s] returned by solver.
        """
        return self._t

    @property
    def y(self) -> _ndarray:
        """
        Saved spatiotemporal solution variables [units].

        Returns
        -------
        y : 2D array
            Solution variables [units] returned by solver.
        """
        return self._y

    @property
    def ydot(self) -> _ndarray:
        """
        Saved solution variable time derivates (i.e., dy/dt) [units].

        Returns
        -------
        ydot : 2D array
            Solution variable time derivatives [units] returned by solver.
        """
        return self._ydot

    @property
    def success(self) -> bool:
        """
        Solver overall exit status.

        Returns
        -------
        success : bool
            ``True`` if no errors, ``False`` otherwise.
        """
        return self._success

    @property
    def onroot(self) -> bool:
        """
        Solver root(event) exit status.

        Returns
        -------
        onroot : bool
            ``True`` if exit on root(event), ``False`` otherwise.
        """
        return self._onroot

    @property
    def message(self) -> str:
        """
        Solver exit message.

        Returns
        -------
        message : str
            Exit message from Sundials IDA solver.
        """
        return self._message

    def solvetime(self, units: str = 's') -> str:
        """
        Print solve time (not including pre/post processing).

        Parameters
        ----------
        units : str, optional
            Units for printed time (``'s'``, ``'min'``, or ``'h'``). The
            default is ``'s'``.

        Returns
        -------
        solvetime : str
            Time for Sundials IDA to solve problem in [units].
        """

        converter = {'s': lambda t: t,
                     'min': lambda t: t / 60.,
                     'h': lambda t: t / 3600.
                     }

        time = converter[units](self._solvetime)

        return f'Solve time: {time:.3f} ' + units

    def ida_fill(self, sol: object, solvetime: float) -> None:
        """
        Fill the instance attributes using the SolverReturn object from the
        Sundials IDA solver.

        Parameters
        ----------
        sol : object
            Sundials IDA SolverReturn object with the following attributes:

            =========== =============================================
            Attribute   Value (type)
            =========== =============================================
            t           solution times (1D array)
            y           solution variables (2D array)
            ydot        solution variable time derivatives (2D array)
            success     overall solver exit status (bool)
            onroot      onroot solver exit status (bool)
            message     solver exit message (str)
            =========== =============================================

        solvetime : float
            Solver integration time, in seconds.

        Returns
        -------
        None.
        """

        self._t = sol.values.t
        self._y = sol.values.y
        self._ydot = sol.values.ydot

        self._success = False
        if sol.flag >= 0:
            self._success = True

        self._onroot = False
        if not isinstance(sol.roots.t, type(None)):
            self._onroot = True

        self._message = sol.message

        self._solvetime = solvetime

    def dict_fill(self, sol: dict) -> None:
        """
        Fill the instance attributes using a dictionary.

        Parameters
        ----------
        sol : dict
            Solution dictionary with the following key/value pairs:

            =========== =============================================
            Key         Value (type)
            =========== =============================================
            t           solution times (1D array)
            y           solution variables (2D array)
            ydot        solution variable time derivatives (2D array)
            success     overall solver exit status (bool)
            onroot      onroot solver exit status (bool)
            message     solver exit message (str)
            solvetime   solver integration time, in seconds (float)
            =========== =============================================

        Returns
        -------
        None.
        """

        self._t = sol['t']
        self._y = sol['y']
        self._ydot = sol['ydot']
        self._success = sol['success']
        self._onroot = sol['onroot']
        self._message = sol['message']
        self._solvetime = sol['solvetime']

    def report(self) -> None:
        """
        Prints the experiment details and solution success report.

        Returns
        -------
        None.
        """

        experiment = 'Experiment('
        for i, (k, v) in enumerate(self._exp.items()):

            if i == 0:
                spacer = ''
            else:
                spacer = ' ' * 11

            if i == len(self._exp) - 1:
                endline = ')'
            else:
                endline = ',\n'

            experiment += spacer + f'{k} = {v}' + endline

        print(experiment + '\n')

        converter = {True: lambda t: None,
                     False: lambda t: str(round(t, 3)) + ' s'
                     }

        solvetime = converter[self._solvetime is None](self._solvetime)

        print(f'Solution(classname = {self.classname},\n'
              f'         success = {self._success},\n'
              f'         onroot = {self._onroot},\n'
              f'         message = {self._message},\n'
              f'         solvetime = {solvetime})'
              )

    def to_dict(self) -> dict:
        """
        Output a dictionary with key/value pairs corresponding to the instance
        attributes and values listed below.

        Returns
        -------
        sol : dict
            Solution dictionary with the following key/value pairs:

            =========== =============================================
            Key         Value (type)
            =========== =============================================
            t           solution times (1D array)
            y           solution variables (2D array)
            ydot        solution variable time derivatives (2D array)
            success     overall solver exit status (bool)
            onroot      onroot solver exit status (bool)
            message     solver exit message (str)
            solvetime   solver integration time, in seconds (float)
            =========== =============================================
        """

        sol = {}
        sol['t'] = self._t
        sol['y'] = self._y
        sol['ydot'] = self._ydot
        sol['success'] = self._success
        sol['onroot'] = self._onroot
        sol['message'] = self._message
        sol['solvetime'] = self._solvetime

        return sol

    def plot(self, *args) -> None:

        if 'current' in args:
            from ..postutils import current
            current(self)

        if 'voltage' in args:
            from ..postutils import voltage
            voltage(self)

        if 'power' in args:
            from ..postutils import power
            power(self)

        if 'IVP' in args:
            from ..postutils import IVP
            IVP(self)

        if 'potentials' in args:
            from ..postutils import potentials
            potentials(self)

        if 'electrolyte' in args:
            from ..postutils import electrolyte
            electrolyte(self)

        if 'intercalation' in args:
            from ..postutils import intercalation
            intercalation(self)

        if 'contours' in args:
            from ..postutils import contours
            contours(self)

    def slice_and_save(self, savename: str, overwrite: bool = False) -> None:
        """
        Save a ``.npz`` file with all spatial, time, and state variables
        separated into 1D, 2D, and 3D arrays. The keys are given below.
        The index order of the 2D and 3D arrays is given with the value
        descriptions.

        ========= ==========================================================
        Key       Value [units] (type)
        ========= ==========================================================
        x_a       x mesh in anode [m] (1D array)
        x_s       x mesh in separator [m] (1D array)
        x_c       x mesh in cathode [m] (1D array)
        x         stacked x mesh for an, sep, and ca [m] (1D array)
        r_a       r mesh for anode particles [m] (1D array)
        r_c       r mesh for cathode particles [m] (1D array)
        t         saved solution times [s] (1D array)
        phie_a    electrolyte potentials at t, x_a [V] (2D array)
        phis_a    electrode potentials at t, x_a [V] (2D array)
        ce_a      electrolyte Li+ at t, x_a [kmol/m^3] (2D array)
        cs_a      electrode Li at t, x_a, r_a [kmol/m^3] (3D array)
        phie_s    electrolyte potentials at t, x_s [V] (2D array)
        ce_s      electrolyte Li+ at t, x_s [kmol/m^3] (2D array)
        phie_c    electrolyte potentials at t, x_c [V] (2D array)
        phis_c    electrode potentials at t, x_c [V] (2D array)
        ce_c      electrolyte Li+ at t, x_c [kmol/m^3] (2D array)
        cs_c      electrode Li at t, x_c, r_c [kmol/m^3] (3D array)
        phie      electrolyte potentials at t, x [V] (2D array)
        ce        electrolyte Li+ at t, x [kmol/m^3] (2D array)
        ie        electrolyte current at t, x boundarys [A/m^2] (2D array)
        j_a       anode Faradaic current at t, x_a [kmol/m^2/s] (2D array)
        j_c       cathode Faradaic current at t, x_c [kmol/m^2/s] (2D array)
        ========= ==========================================================

        Parameters
        ----------
        savename : str
            Either a file name or the absolute/relative file path. The ``.npz``
            extension will be added to the end of the string if it is not
            already there. If only the file name is given, the file will be
            saved in the user's current working directory.

        overwrite : bool, optional
            A flag to overwrite and existing ``.npz`` file with the same name
            if one exists. The default is ``False``.

        Returns
        -------
        None.
        """

        import os

        import numpy as np

        if len(self.postvars) == 0:
            self.post()

        if os.path.exists(savename) and not overwrite:
            raise Exception('save_and_slice file already exists. Overwrite'
                            ' with flag or delete the file and try again.')

        sim = self._sim

        x_a = sim.an.x
        x_s = sim.sep.x
        x_c = sim.ca.x

        x = np.hstack([x_a, x_s, x_c])

        r_a = sim.an.r
        r_c = sim.ca.r

        t = self.t

        phie_a = self.y[:, sim.an.x_ptr('phi_el')]
        phis_a = self.y[:, sim.an.x_ptr('phi_ed')]
        ce_a = self.y[:, sim.an.x_ptr('Li_el')]

        phie_s = self.y[:, sim.sep.x_ptr('phi_el')]
        ce_s = self.y[:, sim.sep.x_ptr('Li_el')]

        phie_c = self.y[:, sim.ca.x_ptr('phi_el')]
        phis_c = self.y[:, sim.ca.x_ptr('phi_ed')]
        ce_c = self.y[:, sim.ca.x_ptr('Li_el')]

        phie = np.hstack([phie_a, phie_s, phie_c])
        ce = np.hstack([ce_a, ce_s, ce_c])

        cs_a = np.zeros([t.size, x_a.size, sim.an.Nr])
        for k in range(sim.an.Nr):
            cs_a[:, :, k] = self.y[:, sim.an.x_ptr('Li_ed', k)] * sim.an.Li_max

        cs_c = np.zeros([t.size, x_c.size, sim.ca.Nr])
        for k in range(sim.ca.Nr):
            cs_c[:, :, k] = self.y[:, sim.ca.x_ptr('Li_ed', k)] * sim.ca.Li_max

        ie = self.postvars['i_el_x']
        j_a = self.postvars['sdot_an']
        j_c = self.postvars['sdot_ca']

        np.savez(savename, x_a=x_a, x_s=x_s, x_c=x_c, x=x, r_a=r_a,
                 r_c=r_c, t=t, phie_a=phie_a, phis_a=phis_a, ce_a=ce_a,
                 phie_s=phie_s, ce_s=ce_s, phie_c=phie_c, phis_c=phis_c,
                 ce_c=ce_c, phie=phie, ce=ce, cs_a=cs_a, cs_c=cs_c,
                 ie=ie, j_a=j_a, j_c=j_c)
