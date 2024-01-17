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

        self._sim = sim
        self._exp = exp

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
        if type(sol.roots.t) != type(None):
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

        spacer = {True: '',
                  False: ' ' * 11
                  }

        endline = {True: ')',
                   False: ',\n'
                   }

        last = len(self._exp)

        experiment = 'Experiment('
        for i, (k, v) in enumerate(self._exp.items()):
            experiment += spacer[i == 0] + f'{k} = {v}' + endline[i == last]

        print(experiment + '\n')

        converter = {True: lambda t: None,
                     False: lambda t: str(round(t, 3)) + ' s'
                     }

        solvetime = converter[self._solvetime == None](self._solvetime)

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

    # def save(self, savename: str, overwrite: bool = False) -> None:
    #     """


    #     Parameters
    #     ----------
    #     savename : str
    #         DESCRIPTION.

    #     overwrite : bool, optional
    #         DESCRIPTION. The default is False.

    #     Raises
    #     ------
    #     Exception
    #         DESCRIPTION.

    #     Returns
    #     -------
    #     None.
    #     """

    #     import os, shutil
    #     from pathlib import Path

    #     import numpy as np
    #     from ruamel.yaml import YAML

    #     savepath = Path(savename)
    #     if os.path.exists(savepath) and not overwrite:
    #         raise Exception(savename + ' already exists. Overwrite with flag.')

    #     elif os.path.exists(savepath):
    #         shutil.rmtree(savepath)

    #     os.mkdir(savepath)

    #     yaml = YAML()
    #     yaml.indent(mapping=4)

    #     path = os.path.dirname(__file__) + '/../default_sims/default_P2D.yaml'

    #     yamlpath = Path(path)
    #     default_yaml = yaml.load(yamlpath)

    #     k_bat = default_yaml['battery'].keys()
    #     k_el = default_yaml['electrolyte'].keys()
    #     k_an = default_yaml['anode'].keys()
    #     k_sep = default_yaml['separator'].keys()
    #     k_ca = default_yaml['cathode'].keys()

    #     sim = self._sim

    #     simdata = {}
    #     simdata['battery'] = dict((k, getattr(sim.bat, k)) for k in k_bat)
    #     simdata['electrolyte'] = dict((k, getattr(sim.el, k)) for k in k_el)
    #     simdata['anode'] = dict((k, getattr(sim.an, k)) for k in k_an)
    #     simdata['separator'] = dict((k, getattr(sim.sep, k)) for k in k_sep)
    #     simdata['cathode'] = dict((k, getattr(sim.ca, k)) for k in k_ca)

    #     yaml.dump(simdata, savepath.joinpath(savename + '_sim.yaml'))

    #     exp = self._exp

    #     yaml.dump(exp, savepath.joinpath(savename + '_exp.yaml'))

    #     soldata = self.to_dict()
    #     soldata['classname'] = self.classname

    #     t = soldata.pop('t')
    #     y = soldata.pop('y')
    #     ydot = soldata.pop('ydot')

    #     yaml.dump(soldata, savepath.joinpath(savename + '_sol.yaml'))

    #     np.savez(savepath.joinpath(savename + '_sol.npz'), t=t, y=y, ydot=ydot)

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
        ie        electrolyte current at t, x boundary [A/m^2] (2D array)
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
            DESCRIPTION. The default is ``False``.

        Returns
        -------
        None.
        """

        import numpy as np

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
            cs_a[:,:,k] = self.y[:, sim.an.x_ptr('Li_ed', k)]*sim.an.Li_max

        cs_c = np.zeros([t.size, x_c.size, sim.ca.Nr])
        for k in range(sim.ca.Nr):
            cs_c[:,:,k] = self.y[:, sim.ca.x_ptr('Li_ed', k)]*sim.ca.Li_max

        ie = self.postvars['i_el_x']
        j_a = self.postvars['sdot_an']
        j_c = self.postvars['sdot_ca']

        np.savez(savename, x_a=x_a, x_s=x_s, x_c=x_c, x=x, r_a=r_a,
                 r_c=r_c, t=t, phie_a=phie_a, phis_a=phis_a, ce_a=ce_a,
                 phie_s=phie_s, ce_s=ce_s, phie_c=phie_c, phis_c=phis_c,
                 ce_c=ce_c, phie=phie, ce=ce, cs_a=cs_a, cs_c=cs_c,
                 ie=ie, j_a=j_a, j_c=j_c)
