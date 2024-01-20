"""
Pseudo-2D Model Package
-----------------------
A packaged pseudo-2D (P2D) model. Build a model using the ``Simulation`` class,
and run any available experiments using its "run" methods, e.g. ``run_CC()``.
The experiments return ``Solution`` class instances with post processing,
plotting, and saving methods.
"""

from . import builder
from . import dae
from . import solutions
from . import postutils


class Simulation(object):

    __slots__ = ['_yamlfile', '_yamlpath', '_flags', 'bat', 'el', 'an', 'sep',
                 'ca', 'sv_0', 'svdot_0', 'lband', 'uband', 'algidx']

    def __init__(self, yamlfile: str = 'default_P2D') -> None:
        """
        Make a P2D simulation capable of running various experiments.

        The initialization will add all of the battery attributes from the
        ``.yaml`` file under its ``bat``, ``el``, ``an``, ``sep``, and ``ca``
        attributes. The ``pre()`` method runs at the end of the initialization
        to add dependent parameters, including the mesh, algebraic indices,
        etc. to the simulation instance. This only happens in ``__init__``,
        which has some implications if the user modifies parameters after
        initialization (see the warning below).

        Parameters
        ----------
        yamlfile : str, optional
            An absolute or relative path to the ``.yaml`` file that defines the
            battery properties. The ``.yaml`` extension will be added to the
            end of the string if it is not already there. The default is
            ``'default_P2D'``, which loads an internal file from the ``bmlite``
            library.

        Warning
        -------
        The user may choose to modify parameters after loading in a ``.yaml``
        file, however, they will need to manually re-run the ``pre()`` method
        if they do so. Otherwise, the dependent parameters may not be
        consistent with the user-defined inputs.

        Notes
        -----
        * The ``'default_P2D.yaml'`` file can be accessed via the
          ``templates()`` method. Please reference this for help building your
          own ``.yaml`` file.
        """

        import os
        from pathlib import Path

        from ruamel.yaml import YAML

        from .builder import Battery, Electrolyte, Electrode, Separator

        if '.yaml' not in yamlfile:
            yamlfile += '.yaml'

        defaults = os.listdir(os.path.dirname(__file__) + '/default_sims')
        if yamlfile in defaults:
            path = os.path.dirname(__file__) + '/default_sims/' + yamlfile
            print('[BatMods WARNING] P2D Simulation: Using a default yaml')
            yamlpath = Path(path)

        elif os.path.exists(yamlfile):
            yamlpath = Path(yamlfile)

        self._yamlfile = yamlfile
        self._yamlpath = yamlpath

        yaml = YAML(typ='safe')
        yamldict = yaml.load(yamlpath)

        self.bat = Battery(**yamldict['battery'])
        self.el = Electrolyte(**yamldict['electrolyte'])
        self.an = Electrode(**yamldict['anode'])
        self.sep = Separator(**yamldict['separator'])
        self.ca = Electrode(**yamldict['cathode'])

        # Function output flags
        self._flags = {}
        self._flags['band'] = False
        self._flags['post'] = False
        self._flags['BC'] = None

        # Pre process dependent parameters, mesh, etc.
        self.pre()

    def pre(self) -> None:
        """
        Pre-process the dependent parameters.

        The dependent parameters include ``A_s``, ``eps_s``, ``eps_AM``,
        ``sigma_s``, and setting the ``material`` classes for each domain. In
        addition, this method determines the mesh, pointers, algebraic indices,
        bandwidth, and initial solution. ``pre()`` is automatically executed
        in the ``__init__()`` method which has some implications if the user
        modifies parameters after initialization (see the warning below).

        Returns
        -------
        None.

        Warning
        -------
        The user may choose to modify parameters after loading in a ``.yaml``
        file, however, they will need to manually re-run the ``pre()`` method
        if they do so. Otherwise, the dependent parameters may not be
        consistent with the user-defined inputs.
        """

        import numpy as np

        from .dae import bandwidth

        # Update dependent parameters
        self.bat.update()
        self.el.update()
        self.an.update()
        self.sep.update()
        self.ca.update()

        # Make meshes
        self.an.make_mesh()

        self.sep.make_mesh()
        self.sep.x = self.an.thick + self.sep.x
        self.sep.xm = self.an.thick + self.sep.xm
        self.sep.xp = self.an.thick + self.sep.xp

        self.ca.make_mesh()
        self.ca.x = self.an.thick + self.sep.thick + self.ca.x
        self.ca.xm = self.an.thick + self.sep.thick + self.ca.xm
        self.ca.xp = self.an.thick + self.sep.thick + self.ca.xp

        # Make pointers
        self.an.make_ptr()

        self.sep.make_ptr()
        for key, val in self.sep.ptr.items():
            if 'shift' not in key and 'off' not in key:
                self.sep.ptr[key] += self.an.ptr['shift']

        self.ca.make_ptr()
        for key, val in self.ca.ptr.items():
            if 'shift' not in key and 'off' not in key:
                self.ca.ptr[key] += self.an.ptr['shift'] + \
                    self.sep.ptr['shift']

        # Initialize potentials [V]
        self.an.phi_0 = 0.

        self.el.phi_0 = -self.an.get_Eeq(self.an.x_0, self.bat.temp)

        self.ca.phi_0 = self.ca.get_Eeq(self.ca.x_0, self.bat.temp) \
                      - self.an.get_Eeq(self.an.x_0, self.bat.temp)

        # Initialize sv and svdot
        self.sv_0 = np.hstack([self.an.sv_0(self.el), self.sep.sv_0(self.el),
                               self.ca.sv_0(self.el)])

        self.svdot_0 = np.zeros_like(self.sv_0)

        # Algebraic indices
        self.algidx = list(self.an.algidx()) + list(self.sep.algidx()) \
                    + list(self.ca.algidx())

        # Determine the bandwidth
        self.lband, self.uband, _ = bandwidth(self)

    def j_pattern(self) -> None:
        """
        Plot the Jacobian pattern.

        Runs the ``bmlite.P2D.dae.bandwidth`` function to determine and plot
        the Jacobian pattern.

        Returns
        -------
        None.

        See also
        --------
        bmlite.P2D.dae.bandwidth
        """

        import matplotlib.pyplot as plt

        from ..plotutils import format_ticks, show
        from .dae import bandwidth

        lband, uband, j_pat = bandwidth(self)

        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=[5, 3.5])

        ax.spy(j_pat)
        ax.text(0.1, 0.2, 'lband: ' + str(lband), transform=ax.transAxes)
        ax.text(0.1, 0.1, 'uband: ' + str(uband), transform=ax.transAxes)

        format_ticks(ax)
        show(fig)

    def run_CC(self, exp: dict, **kwargs) -> object:
        """
        Runs a constant current experiment specified by the details given in
        the experiment dictionary ``exp``.

        Parameters
        ----------
        exp : dict
            The constant current experimental details. Required keys and
            descriptions are listed below:

            =========== ==============================================
            Key         Value [units] (type)
            =========== ==============================================
            C_rate      C-rate (+ charge, - discharge) [1/h] (*float*)
            t_min       minimum time [s] (*float*)
            t_max       maximum time [s] (*float*)
            Nt          number of time discretizations [-] (*int*)
            =========== ==============================================

        **kwargs : dict, optional
            The keyword arguments specify the Sundials IDA solver options. A
            partial list of options/defaults is given below:

            =============== =================================================
            Key             Description (type or options, default)
            =============== =================================================
            rtol            relative tolerance (*float*, 1e-6)
            atol            absolute tolerance (*float*, 1e-9)
            linsolver       linear solver (``{'dense', 'band'}``, ``'band'``)
            lband           width of the lower band (*int*, ``self.lband``)
            uband           width of the upper band (*int*, ``self.uband``)
            max_step_size   maximum time step (*float*, 0. -> unrestricted)
            rootfn          root/event function (*Callable*, ``None``)
            nr_rootfns      number of events in ``'rootfn'`` (*int*, 0)
            =============== =================================================

        Returns
        -------
        sol : CCSolution object
            Solution class with the returned variable values, messages, exit
            flags, etc. from the IDA solver. The returned ``CCSolution``
            instance includes post processing, plotting, and saving methods.

        See also
        --------
        bmlite.IDASolver
        bmlite.P2D.solutions.CCSolution
        """

        import time

        import numpy as np

        from .. import IDASolver
        from .dae import residuals
        from .solutions import CCSolution

        # Boundary condition flag
        self._flags['BC'] = 'current'

        # Initialize solution
        sol = CCSolution(self, exp)

        # Solver options
        options = {'user_data': (self, exp),
                   'linsolver': 'band',
                   'lband': self.lband,
                   'uband': self.uband,
                   'algidx': self.algidx
                   }

        for key, val in kwargs.items():
            options[key] = val

        solver = IDASolver(residuals, **options)

        # Run the simulation
        t_span = np.linspace(exp['t_min'], exp['t_max'], exp['Nt'])

        start = time.time()
        idasol = solver.solve(t_span, self.sv_0, self.svdot_0)
        solvetime = time.time() - start

        sol.ida_fill(idasol, solvetime)

        return sol

    def run_CV(self, exp: dict, **kwargs) -> object:
        """
        Runs a constant voltage experiment specified by the details given in
        the experiment dictionary ``exp``.

        Parameters
        ----------
        exp : dict
            The constant voltage experimental details. Required keys and
            descriptions are listed below:

            =========== ==========================================
            Key         Value [units] (type)
            =========== ==========================================
            V_ext       externally applied voltage [V] (*float*)
            t_min       minimum time [s] (*float*)
            t_max       maximum time [s] (*float*)
            Nt          number of time discretizations [-] (*int*)
            =========== ==========================================

        **kwargs : dict, optional
            The keyword arguments specify the Sundials IDA solver options. A
            partial list of options/defaults is given below:

            =============== =================================================
            Key             Description (type or options, default)
            =============== =================================================
            rtol            relative tolerance (*float*, 1e-6)
            atol            absolute tolerance (*float*, 1e-9)
            linsolver       linear solver (``{'dense', 'band'}``, ``'band'``)
            lband           width of the lower band (*int*, ``self.lband``)
            uband           width of the upper band (*int*, ``self.uband``)
            max_step_size   maximum time step (*float*, 0. -> unrestricted)
            rootfn          root/event function (*Callable*, ``None``)
            nr_rootfns      number of events in ``'rootfn'`` (*int*, 0)
            =============== =================================================

        Returns
        -------
        sol : CVSolution object
            Solution class with the returned variable values, messages, exit
            flags, etc. from the IDA solver. The returned ``CVSolution``
            instance includes post processing, plotting, and saving methods.

        See also
        --------
        bmlite.IDASolver
        bmlite.P2D.solutions.CVSolution
        """

        import time

        import numpy as np

        from .. import IDASolver
        from .dae import residuals
        from .solutions import CVSolution

        # Boundary condition flag
        self._flags['BC'] = 'voltage'

        # Initialize solution
        sol = CVSolution(self, exp)

        # Solver options
        options = {'user_data': (self, exp),
                   'linsolver': 'band',
                   'lband': self.lband,
                   'uband': self.uband,
                   'algidx': self.algidx
                   }

        for key, val in kwargs.items():
            options[key] = val

        solver = IDASolver(residuals, **options)

        # Run the simulation
        t_span = np.linspace(exp['t_min'], exp['t_max'], exp['Nt'])

        start = time.time()
        idasol = solver.solve(t_span, self.sv_0, self.svdot_0)
        solvetime = time.time() - start

        sol.ida_fill(idasol, solvetime)

        if not sol.success:
            print('[BatMods WARNING] run_CV: bad initstep, trying to resolve')

            V_ext = exp['V_ext']
            V_init = self.sv_0[self.ca.ptr['phi_ed']]

            start = time.time()

            for wt in [0.9, 0.8, 0.7, 0.5, 0.4, 0.3, 0.2, 0.1]:

                exp['V_ext'] = wt * V_init + (1 - wt) * V_ext
                init = solver.init_step(exp['t_min'], self.sv_0, self.svdot)

                exp['V_ext'] = V_ext
                idasol = solver.solve(t_span, init.values.y, init.values.ydot)

                if idasol.flag >= 0:
                    break

            if not idasol.flag >= 0:
                print('[BatMods ERROR] run_CV: failed to resolve bad initstep')
            else:
                print('[BatMods NOTE] run_CV: initstep successfully resolved')

            solvetime = time.time() - start

            sol.ida_fill(idasol, solvetime)

        # Reset boundary condition flag
        self._flags['BC'] = None

        return sol

    def run_CP(self, exp: dict, **kwargs) -> object:
        """
        Runs a constant power experiment specified by the details given in
        the experiment dictionary ``exp``.

        Parameters
        ----------
        exp : dict
            The constant power experimental details. Required keys and
            descriptions are listed below:

            ======== ========================================================
            Key      Value [units] (type)
            ======== ========================================================
            P_ext    external power (+ charge, - discharge) [W/m^2] (*float*)
            t_min    minimum time [s] (*float*)
            t_max    maximum time [s] (*float*)
            Nt       number of time discretizations [-] (*int*)
            ======== ========================================================

        **kwargs : dict, optional
            The keyword arguments specify the Sundials IDA solver options. A
            partial list of options/defaults is given below:

            =============== =================================================
            Key             Description (type or options, default)
            =============== =================================================
            rtol            relative tolerance (*float*, 1e-6)
            atol            absolute tolerance (*float*, 1e-9)
            linsolver       linear solver (``{'dense', 'band'}``, ``'band'``)
            lband           width of the lower band (*int*, ``self.lband``)
            uband           width of the upper band (*int*, ``self.uband``)
            max_step_size   maximum time step (*float*, 0. -> unrestricted)
            rootfn          root/event function (*Callable*, ``None``)
            nr_rootfns      number of events in ``'rootfn'`` (*int*, 0)
            =============== =================================================

        Returns
        -------
        sol : CPSolution object
            Solution class with the returned variable values, messages, exit
            flags, etc. from the IDA solver. The returned ``CPSolution``
            instance includes post processing, plotting, and saving methods.

        See also
        --------
        bmlite.IDASolver
        bmlite.P2D.solutions.CPSolution
        """

        import time

        import numpy as np

        from .. import IDASolver
        from .dae import residuals
        from .solutions import CPSolution

        # Boundary condition flag
        self._flags['BC'] = 'power'

        # Initialize solution
        sol = CPSolution(self, exp)

        # Solver options
        options = {'user_data': (self, exp),
                   'linsolver': 'band',
                   'lband': self.lband,
                   'uband': self.uband,
                   'algidx': self.algidx
                   }

        for key, val in kwargs.items():
            options[key] = val

        solver = IDASolver(residuals, **options)

        # Run the simulation
        t_span = np.linspace(exp['t_min'], exp['t_max'], exp['Nt'])

        start = time.time()
        idasol = solver.solve(t_span, self.sv_0, self.svdot_0)
        solvetime = time.time() - start

        sol.ida_fill(idasol, solvetime)

        # Reset boundary condition flag
        self._flags['BC'] = None

        return sol

    def copy(self) -> object:
        from copy import deepcopy

        return deepcopy(self)


def templates(sim: str | int = None, exp: str | int = None) -> None:
    """
    Print simulation and/or experiment templates. If both ``sim`` and ``exp``
    are ``None``, a list of available templates will be printed. Otherwise, if
    a name or index is given, that template will print to the console.

    Parameters
    ----------
    sim : str | int, optional
        Simulation template file name or index. The default is ``None``.

    exp : str | int, optional
        Experiment template file name or index. The default is ``None``.

    Returns
    -------
    None.
    """

    from .. import _templates

    _templates(__file__, 'P2D', sim, exp)
