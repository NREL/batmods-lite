from numpy import ndarray as _ndarray


class BaseSolution(object):
    """
    Base methods for all P2D Solution classes.
    """

    __slots__ = ['_sim', '_exp', '_t', '_y', '_ydot', '_success', '_onroot',
                 '_message', '_solvetime', 'postvars']

    def __init__(self, sim: object, exp: dict) -> None:
        """
        When initialized, a copy of the SPM simulation object and experimental
        details for the solution instance are stored, all other instance
        attributes are set to ``None``. These are filled later with a "fill"
        method (e.g., ``ida_fill()`` or ``dict_fill()``).

        Parameters
        ----------
        sim : SPM Simulation object
            The SPM Simulation instance used to produce the solution.

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

        self.postvars = {}

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
            ``True`` if exit on root/event, ``False`` otherwise.
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
            Units for printed time ``{'s', 'min', 'h'}``. The default is
            ``'s'``.

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

            =========== ========================================
            Attribute   Value [units] (*type*)
            =========== ========================================
            t           solution times [s] (*1D array*)
            y           solution variables [units] (*2D array*)
            ydot        dy/dt derivatives [units/s] (*2D array*)
            success     overall solver exit status (*bool*)
            onroot      onroot solver exit status (*bool*)
            message     solver exit message (*str*)
            =========== ========================================

        solvetime : float
            Solver integration time, in seconds.

        Returns
        -------
        None.
        """

        self._t = sol.values.t
        self._y = sol.values.y
        self._ydot = sol.values.ydot
        self._success = bool(sol.flag >= 0)
        self._onroot = bool(not isinstance(sol.roots.t, type(None)))
        self._message = sol.message
        self._solvetime = solvetime

    def dict_fill(self, sol: dict) -> None:
        """
        Fill the instance attributes using a dictionary.

        Parameters
        ----------
        sol : dict
            Solution dictionary with the following key/value pairs:

            =========== ========================================
            Key         Value [units] (*type*)
            =========== ========================================
            t           solution times [s] (*1D array*)
            y           solution variables [units] (*2D array*)
            ydot        dy/dt derivatives [units/s] (*2D array*)
            success     overall solver exit status (*bool*)
            onroot      onroot solver exit status (*bool*)
            message     solver exit message (*str*)
            solvetime   solver integration time [s] (*float*)
            =========== ========================================

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
              f'         solvetime = {solvetime})\n'
              )

    def _save_dict(self) -> dict:
        """
        Output a dictionary with key/value pairs corresponding to the instance
        attributes and values listed below.

        Returns
        -------
        sol : dict
            Solution dictionary with the following key/value pairs:

            =========== ========================================
            Key         Value [units] (*type*)
            =========== ========================================
            t           solution times [s] (*1D array*)
            y           solution variables [units] (*2D array*)
            ydot        dy/dt derivatives [units/s] (*2D array*)
            success     overall solver exit status (*bool*)
            onroot      onroot solver exit status (*bool*)
            message     solver exit message (*str*)
            solvetime   solver integration time [s] (*float*)
            =========== ========================================
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

    def post(self) -> None:
        from ..postutils import post

        self.postvars = post(self)

    def plot(self, *args: str) -> None:
        """
        Generates requested plots based on ``*args``.

        Parameters
        ----------
        *args : str
            Use any number of the following arguments to see the described
            plots:

            ================= ===============================================
            arg               Plot description
            ================= ===============================================
            'current'         current density [A/m^2] vs. time [s]
            'voltage'         cell voltage [V] vs. time [s]
            'power'           power density [W/m^2] vs. time [s]
            'ivp'             combined current, voltage, and power
            'potentials'      anode, cathode, and electrolyte potentials [V]
            'intercalation'   anode/cathode particle Li fractions vs. r and t
            'pixels'          pixel plots for most 2D variables
            ================= ===============================================

        Returns
        -------
        None.
        """

        if len(self.postvars) == 0:
            self.post()

        if 'current' in args:
            from ..postutils import current
            current(self)

        if 'voltage' in args:
            from ..postutils import voltage
            voltage(self)

        if 'power' in args:
            from ..postutils import power
            power(self)

        if 'ivp' in args:
            from ..postutils import ivp
            ivp(self)

        if 'potentials' in args:
            from ..postutils import potentials
            potentials(self)

        if 'intercalation' in args:
            from ..postutils import intercalation
            intercalation(self)

        if 'pixels' in args:
            from ..postutils import pixels
            pixels(self)

    def to_dict(self) -> dict:
        """
        Creates a dict with all spatial, time, and state variables separated
        into 1D and 2D arrays. The keys are given below.

        ========= =======================================================
        Key       Value [units] (*type*)
        ========= =======================================================
        r_a       r mesh for anode particles [m] (*1D array*)
        r_c       r mesh for cathode particles [m] (*1D array*)
        t         saved solution times [s] (*1D array*)
        phis_a    anode electrode potentials at t [V] (*1D array*)
        cs_a      electrode Li at t, r_a [kmol/m^3] (*2D array*)
        phis_c    cathode electrode potentials at t [V] (*1D array*)
        cs_c      electrode Li at t, r_c [kmol/m^3] (*2D array*)
        phie      electrolyte potentials at t [V] (*1D array*)
        j_a       anode Faradaic current at t [kmol/m^2/s] (*1D array*)
        j_c       cathode Faradaic current at t [kmol/m^2/s] (*1D array*)
        ========= =======================================================

        Parameters
        ----------
        None.

        Returns
        -------
        sol_dict : dict
            A dictionary containing the solution
        """

        if len(self.postvars) == 0:
            self.post()

        sim = self._sim

        an_sol = sim.an.to_dict(self)
        ca_sol = sim.ca.to_dict(self)
        el_sol = sim.el.to_dict(self)

        sol_dict = {}

        sol_dict['r_a'] = sim.an.r
        sol_dict['r_c'] = sim.ca.r

        sol_dict['t'] = self.t

        sol_dict['phis_a'] = an_sol['phis']
        sol_dict['phis_c'] = ca_sol['phis']
        sol_dict['phie'] = el_sol['phie']

        sol_dict['cs_a'] = an_sol['cs']
        sol_dict['cs_c'] = ca_sol['cs']

        sol_dict['j_a'] = self.postvars['sdot_an']
        sol_dict['j_c'] = self.postvars['sdot_ca']

        return sol_dict

    def save_sliced(self, savename: str, overwrite: bool = False) -> None:
        """
        Save a ``.npz`` file with all spatial, time, and state variables
        separated into 1D and 2D arrays. The keys are given below. The index
        order of the 2D arrays is given with the value descriptions.

        ========= =======================================================
        Key       Value [units] (*type*)
        ========= =======================================================
        r_a       r mesh for anode particles [m] (*1D array*)
        r_c       r mesh for cathode particles [m] (*1D array*)
        t         saved solution times [s] (*1D array*)
        phis_a    anode electrode potentials at t [V] (*1D array*)
        cs_a      electrode Li at t, r_a [kmol/m^3] (*2D array*)
        phis_c    cathode electrode potentials at t [V] (*1D array*)
        cs_c      electrode Li at t, r_c [kmol/m^3] (*2D array*)
        phie      electrolyte potentials at t [V] (*1D array*)
        j_a       anode Faradaic current at t [kmol/m^2/s] (*1D array*)
        j_c       cathode Faradaic current at t [kmol/m^2/s] (*1D array*)
        ========= =======================================================

        Parameters
        ----------
        savename : str
            Either a file name or the absolute/relative file path. The ``.npz``
            extension will be added to the end of the string if it is not
            already there. If only the file name is given, the file will be
            saved in the user's current working directory.

        overwrite : bool, optional
            A flag to overwrite an existing ``.npz`` file with the same name
            if one exists. The default is ``False``.

        Returns
        -------
        None.
        """

        import os

        import numpy as np

        if '.npz' not in savename:
            savename += '.npz'

        if os.path.exists(savename) and not overwrite:
            raise FileExistsError(savename + ' already exists. Use overwrite'
                                  ' flag or delete the file and try again.')

        sol_dict = self.to_dict()

        np.savez(savename, **sol_dict)
