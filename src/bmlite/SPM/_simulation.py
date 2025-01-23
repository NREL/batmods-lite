from __future__ import annotations
from typing import TYPE_CHECKING

import os
import time
from pathlib import Path
from copy import deepcopy

import numpy as np
from ruamel.yaml import YAML
import matplotlib.pyplot as plt

if TYPE_CHECKING:  # pragma: no cover
    from bmlite import Experiment
    from ._solutions import StepSolution, CycleSolution


class Simulation:

    __slots__ = ['_yamlfile', '_yamlpath', '_t0', '_sv0', '_svdot0', '_lband',
                 '_uband', '_algidx', 'c', 'bat', 'el', 'an', 'ca']

    def __init__(self, yamlfile: str = 'graphite_nmc532') -> None:
        """
        Make a SPM simulation capable of running various experiments.

        The initialization will add all of the battery attributes from the
        ``.yaml`` file under its ``bat``, ``el``, ``an``, and ``ca``
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
            ``'default_SPM'``, which loads an internal file from the ``bmlite``
            library.

        Warning
        -------
        The user may choose to modify parameters after loading in a ``.yaml``
        file, however, they will need to manually re-run the ``pre()`` method
        if they do so. Otherwise, the dependent parameters may not be
        consistent with the user-defined inputs.

        See also
        --------
        bmlite.templates :
            Get help making your own ``.yaml`` file by starting with the
            default template.

        """

        from .. import Constants
        from .._utils import short_warn
        from .domains import Battery, Electrolyte, Electrode

        if '.yaml' not in yamlfile:
            yamlfile += '.yaml'

        defaults = os.listdir(os.path.dirname(__file__) + '/templates')
        if yamlfile in defaults:
            path = os.path.dirname(__file__) + '/templates/' + yamlfile
            short_warn(f"SPM Simulation: Using default {yamlfile}")
            yamlpath = Path(path)

        elif os.path.exists(yamlfile):
            yamlpath = Path(yamlfile)

        else:
            raise FileNotFoundError(yamlfile)

        self._yamlfile = yamlfile
        self._yamlpath = yamlpath

        yaml = YAML(typ='safe')
        yamldict = yaml.load(yamlpath)

        self.c = Constants()
        self.bat = Battery(**yamldict['battery'])
        self.el = Electrolyte(**yamldict['electrolyte'])
        self.an = Electrode('anode', **yamldict['anode'])
        self.ca = Electrode('cathode', **yamldict['cathode'])

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

        # Update dependent parameters
        self.bat.update()
        self.el.update()
        self.an.update()
        self.ca.update()

        # Make meshes/pointers
        self.an.make_mesh()
        self.el.make_mesh(pshift=self.an.ptr['shift'])
        self.ca.make_mesh(pshift=self.an.ptr['shift'] + self.el.ptr['shift'])

        # Initialize potentials [V]
        self.an.phi_0 = 0.

        self.el.phi_0 = -self.an.get_Eeq(self.an.x_0)

        self.ca.phi_0 = self.ca.get_Eeq(self.ca.x_0) \
                      - self.an.get_Eeq(self.an.x_0)

        # Initialize sv and svdot
        self._t0 = 0.
        self._sv0 = np.hstack([self.an.sv0(), self.el.sv0(), self.ca.sv0()])
        self._svdot0 = np.zeros_like(self._sv0)

        # Algebraic indices
        self._algidx = self.an.algidx().tolist() + self.el.algidx().tolist() \
                     + self.ca.algidx().tolist()

        # Determine the bandwidth
        # self._lband, self._uband, _ = bandwidth(self)
        self._lband = int(self.el.ptr['phi_el'] - self.an.r_ptr['Li_ed'][-1])
        self._uband = int(self.el.ptr['phi_el'] - self.an.r_ptr['Li_ed'][-1])

    def j_pattern(self, plot: bool = True,
                  return_bands: bool = False) -> tuple[int] | None:
        """
        Determine the Jacobian pattern.

        Parameters
        ----------
        plot : bool, optional
            Whether or not to plot the Jacobian pattern. The default is True.
        return_barnds : bool, optional
            Whether or not to return the half bandwidths (lower, upper). The
            default is False.

        Returns
        -------
        lband : int
            The lower half bandwidth. Only returned if `return_bands=True`.
        uband : int
            The upper half bandwidth. Only returned if `return_bands=True`.

        """

        from .dae import residuals
        from .._utils import ExitHandler
        from .._core._idasolver import bandwidth
        from bmlite.plotutils import format_ticks

        t0 = 0.
        y0 = np.hstack([self.an.sv0(), self.el.sv0(), self.ca.sv0()])
        yp0 = np.zeros_like(y0)

        step = {
            'mode': 'current',
            'units': 'C',
            'value': lambda t: 0.,
        }

        userdata = (self, step)

        lband, uband, j_pat = bandwidth(residuals, t0, y0, yp0, userdata,
                                        return_pattern=True)

        if plot:
            _, ax = plt.subplots(nrows=1, ncols=1, figsize=[4, 4],
                                 layout='constrained')

            ax.spy(j_pat)
            ax.text(0.1, 0.2, 'lband: ' + str(lband), transform=ax.transAxes)
            ax.text(0.1, 0.1, 'uband: ' + str(uband), transform=ax.transAxes)

            format_ticks(ax)

            if not plt.isinteractive():
                ExitHandler.register_atexit(plt.show)

        if return_bands:
            return lband, uband

    def run_step(self, expr: Experiment, stepidx: int) -> StepSolution:
        """
        Run a single experimental step.

        Parameters
        ----------
        expr : Experiment
            An experiment instance.
        stepidx : int
            Step index to run. The first step has index 0.

        Returns
        -------
        :class:`~bmlite.SPM.StepSolution`
            Solution to the experimental step.

        Warning
        -------
        The model's internal state is changed at the end of each experimental
        step. Consequently, you should not run steps out of order. You should
        always start with ``stepidx = 0`` and then progress to the subsequent
        steps afterward. Run ``pre()`` after your last step to reset the state
        back to a rested condition at 'soc0', if needed. Alternatively, you
        can continue running experiments back-to-back without a pre-processing
        in between if you want the following experiment to pick up from the
        same state that the last experiment ended.

        See also
        --------
        Experiment : Build an experiment.
        StepSolution : Wrapper for a single-step solution.

        Notes
        -----
        Using the ``run()`` loops through all steps in an experiment and then
        stitches their solutions together. Most of the time, this is more
        convenient. However, advantages for running step-by-step is that it
        makes it easier to fine tune solver options, and allows for analyses
        or control decisions in the middle of an experiment.

        """

        from bmlite import IDASolver

        from .dae import residuals
        from ._solutions import StepSolution

        step = expr.steps[stepidx].copy()
        options = expr._step_options[stepidx].copy()

        if not callable(step['value']):
            value = step['value']
            step['value'] = lambda t: value

        options['userdata'] = (self, step)
        options['calc_initcond'] = 'yp0'
        options['algebraic_idx'] = self._algidx
        options['linsolver'] = 'band'
        options['lband'] = self._lband
        options['uband'] = self._uband

        if step['limits'] is not None:
            _setup_eventsfn(step['limits'], options)

        solver = IDASolver(residuals, **options)

        start = time.time()
        idasoln = solver.solve(step['tspan'], self._sv0, self._svdot0)
        timer = time.time() - start

        soln = StepSolution(self, idasoln, timer)

        self._t0 = soln.t[-1]
        self._sv0 = soln.y[-1].copy()
        self._svdot0 = soln.yp[-1].copy()

        return soln

    def run(self, expr: Experiment, reset_state: bool = True,
            t_shift: float = 1e-3, bar: bool = False) -> CycleSolution:
        """
        Run a full experiment.

        Parameters
        ----------
        expr : Experiment
            An experiment instance.
        reset_state : bool, optional
            If True (default), the internal state of the model will be reset
            back to a rested condition at 'soc0' at the end of all steps. When
            False, the state does not reset. Instead it will update to match
            the final state of the last experimental step.
        t_shift : float, optional
            Time (in seconds) to shift step solutions by when stitching them
            together. If zero the end time of each step overlaps the starting
            time of its following step. The default is 1e-3.
        bar : bool, optional
            Displays a progress bar showing percentage of completed steps when
            True. The default is False.

        Returns
        -------
        :class:`~bmlite.SPM.CycleSolution`
            A stitched solution with all experimental steps.

        Warning
        -------
        The default behavior resets the model's internal state back to a rested
        condition at 'soc0' by calling the ``pre()`` method at the end of all
        steps. This means that if you run a second experiment afterward, it
        will not start where the previous one left off. Instead, it will start
        from the original rested condition that the model initialized with. You
        can bypass this by using ``reset_state=False``, which keeps the state
        at the end of the final experimental step.

        See also
        --------
        Experiment : Build an experiment.
        CycleSolution : Wrapper for an all-steps solution.

        """

        from .._utils import ProgressBar
        from ._solutions import CycleSolution

        progbar = ProgressBar(expr.num_steps)

        solns = []
        for i in range(expr.num_steps):
            solns.append(self.run_step(expr, i))

            if bar:
                progbar.update(i+1)

        soln = CycleSolution(*solns, t_shift=t_shift)

        self._t0 = 0.
        if reset_state:
            self.pre()

        return soln

    def copy(self) -> object:
        """
        Create a copy of the Simulation instance.

        Returns
        -------
        sim : SPM Simulation object
            A unique copy (stored separately in memory) of the Simulation
            instance.

        """
        return deepcopy(self)


class _EventsFunction:
    """Events function callable."""

    def __init__(self, limits: tuple[str, float]) -> None:
        """
        This class is a generalized events function callable.

        Parameters
        ----------
        limits : tuple[str, float]
            A tuple of event function criteria arranged as ordered pairs of
            limit names and values, e.g., ('time_h', 10., 'voltage_V', 4.2).

        """

        self.keys = limits[0::2]
        self.values = limits[1::2]
        self.size = len(self.keys)

    def __call__(self, t: float, sv: np.ndarray, svdot: np.ndarray,
                 events: np.ndarray, inputs: dict) -> None:
        """
        Solver-structured events function.

        The IDASolver requires an events function in this exact form. Rather
        than outputting the events array, the function returns None, but
        fills the 'events' input array with root functions. If any 'events'
        index equals zero, the solver will exit prior to 'tspan[-1]'.

        Parameters
        ----------
        t : float
            Value of time [s].
        sv : 1D np.array
            State variables at time t.
        svdot : 1D np.array
            State variable time derivatives at time t.
        events : 1D np.array
            An array of root/event functions. During integration, the solver
            will exit prior to 'tspan[-1]' if any 'events' index equals zero.
        inputs : dict
            Dictionary detailing an experimental step, with the 'roots' key
            added and filled within the `rhs_funcs()' method.

        Returns
        -------
        None.

        """

        inputs = inputs[1]

        for i, (key, value) in enumerate(zip(self.keys, self.values)):
            events[i] = inputs['events'][key] - value


def _setup_eventsfn(limits: tuple[str, float], kwargs: dict) -> None:
    """
    Set up an events function for the IDASolver.

    The IDASolver requires two keyword arguments to be set when using event
    functions. The first is 'eventsfn' which requires a Callable. The second
    is 'num_events' with allocates memory to an array that stores the events
    function values.

    Parameters
    ----------
    limits : tuple[str, float]
        A tuple of event function criteria arranged as ordered pairs of limit
        names and values, e.g., ('time_h', 10., 'voltage_V', 4.2).
    kwargs : dict
        The IDASolver keyword argumnents dictionary. Both the 'eventsfn' and
        'num_events' keyword arguments must be added to 'kwargs'.

    Returns
    -------
    None.

    """

    eventsfn = _EventsFunction(limits)

    kwargs['eventsfn'] = eventsfn
    kwargs['num_events'] = eventsfn.size
