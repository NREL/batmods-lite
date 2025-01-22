from __future__ import annotations
from typing import Iterable, TYPE_CHECKING

import textwrap
from copy import deepcopy

import atexit
import numpy as np
import matplotlib.pyplot as plt

from bmlite import IDAResult

if TYPE_CHECKING:  # pragma: no cover
    from ._simulation import Simulation


class BaseSolution(IDAResult):
    """Base SPM solution."""

    def __init__(self) -> None:
        """
        The base solution class is a parent class to both the StepSolution
        and CycleSolution classes. Inheriting from this class gives each
        solution instance a 'vars' dictionary, access to the 'plot' method,
        and ensures that the slicing of the solution vector into 'vars' is
        consistent between all solutions.

        Returns
        -------
        None.

        """

        self.vars = {}
        self._postvars = False

    def __repr__(self) -> str:  # pragma: no cover
        """
        Return a readable repr string.

        Returns
        -------
        readable : str
            A console-readable instance representation.

        """

        classname = self.__class__.__name__

        def wrap_string(label: str, value: list, width: int):
            if isinstance(value, Iterable):
                value = list(value)
            else:
                value = [value]

            indent = ' '*(len(label) + 1)

            if classname == 'StepSolution' and len(value) == 1:
                text = label + f"{value[0]!r}"
            else:
                text = label + "[" + ", ".join(f"{v!r}" for v in value) + "]"

            return textwrap.fill(text, width=width, subsequent_indent=indent)

        data = [
            wrap_string('    success=', self.success, 79),
            wrap_string('    status=', self.status, 79),
            wrap_string('    nfev=', self.nfev, 79),
            wrap_string('    njev=', self.njev, 79),
            wrap_string('    vars=', self.vars.keys(), 79),
        ]

        summary = f"    solvetime={self.solvetime},"
        for d in data:
            summary += f"\n{d},"

        readable = f"{classname}(\n{summary}\n)"

        return readable

    def post(self) -> None:
        from .postutils import post
        postvars = post(self)

        self.vars['an']['sdot'] = postvars['sdot_an']
        self.vars['ca']['sdot'] = postvars['sdot_ca']

        self._postvars = True

    def simple_plot(self, x: str, y: str, **kwargs) -> None:
        """
        Plot any two variables in 'vars' against each other.

        Parameters
        ----------
        x : str
            A variable key in 'vars' to be used for the x-axis.
        y : str
            A variable key in 'vars' to be used for the y-axis.
        **kwargs : dict, optional
            Keyword arguments to pass through to `plt.plot()`. For more info
            please refer to documentation for `maplotlib.pyplot.plot()`.

        Returns
        -------
        None.

        """

        plt.figure()
        plt.plot(self.vars[x], self.vars[y], **kwargs)

        variable, units = x.split('_')
        xlabel = variable.capitalize() + ' [' + units + ']'

        variable, units = y.split('_')
        ylabel = variable.capitalize() + ' [' + units + ']'

        plt.xlabel(xlabel)
        plt.ylabel(ylabel)

        if not plt.isinteractive():
            atexit.register(plt.show)

    def complex_plot(self, *args: str) -> None:
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
            'potentials'      anode, cathode, and electrolyte potentials [V]
            'intercalation'   anode/cathode particle Li fractions vs. r and t
            'pixels'          pixel plots for solid Li concentrations
            ================= ===============================================

        Returns
        -------
        None.

        """

        if not self._postvars:
            self.post()

        if 'potentials' in args:
            from .postutils import potentials
            potentials(self)

        if 'intercalation' in args:
            from .postutils import intercalation
            intercalation(self)

        if 'pixels' in args:
            from .postutils import pixels
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

        if not self._postvars:
            self.post()

        vars = {
            'r_a': self.vars['an']['r'],
            'r_c': self.vars['ca']['r'],
            't': self.vars['time_s'],
            'phis_a': self.vars['an']['phis'],
            'cs_a': self.vars['an']['cs'],
            'phis_c': self.vars['ca']['phis'],
            'cs_c': self.vars['ca']['cs'],
            'phie': self.vars['el']['phie'],
            'j_a': self.vars['an']['sdot'],
            'j_c': self.vars['ca']['sdot'],
        }

        return vars

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

    def _fill_vars(self) -> None:
        """
        Fills the 'vars' dictionary by slicing the SolverReturn solution
        states. Users should generally only access the solution via 'vars'
        since names are more intuitive than interpreting 'y' directly.

        Returns
        -------
        None.

        """

        sim = self._sim
        an = sim.an.to_dict(self)
        ca = sim.ca.to_dict(self)
        el = sim.el.to_dict(self)

        # domain variables
        self.vars['an'] = an
        self.vars['ca'] = ca
        self.vars['el'] = el

        # stored time
        time = self.t

        self.vars['time_s'] = time
        self.vars['time_min'] = time / 60.
        self.vars['time_h'] = time / 3600.

        # common variables
        voltage = ca['phis']
        current = sim.ca._boundary_current(self)

        self.vars['current_A'] = current
        self.vars['current_C'] = current / sim.bat.cap
        self.vars['voltage_V'] = voltage
        self.vars['power_W'] = current*voltage


class StepSolution(BaseSolution):
    """Single-step solution."""

    def __init__(self, sim: Simulation, idasoln: IDAResult,
                 timer: float) -> None:
        """
        A solution instance for a single experimental step.

        Parameters
        ----------
        sim : Simulation
            The simulation instance that was run to produce the solution.
        idasoln : IDAResult
            The unformatted solution returned by IDASolver.
        timer : float
            Amount of time it took for IDASolver to perform the integration.

        """

        super().__init__()

        self._sim = sim.copy()

        self.message = idasoln.message
        self.success = idasoln.success
        self.status = idasoln.status

        self.t = idasoln.t
        self.y = idasoln.y
        self.yp = idasoln.yp

        self.i_events = idasoln.i_events
        self.t_events = idasoln.t_events
        self.y_events = idasoln.y_events
        self.yp_events = idasoln.yp_events

        self.nfev = idasoln.nfev
        self.njev = idasoln.njev

        self._timer = timer

        self._fill_vars()

    @property
    def solvetime(self) -> str:
        """
        Print a statement specifying how long IDASolver spent integrating.

        Returns
        -------
        solvetime : str
            An f-string with the solver integration time in seconds.

        """
        return f"{self._timer:.3f} s"


class CycleSolution(BaseSolution):
    """All-step solution."""

    def __init__(self, *soln: StepSolution, t_shift: float = 1e-3) -> None:
        """
        A solution instance with all experiment steps stitch together into
        a single cycle.

        Parameters
        ----------
        *soln : StepSolution
            All unpacked StepSolution instances to stitch together. The given
            steps should be given in the same sequential order that they were
            run.
        t_shift : float
            Time (in seconds) to shift step solutions by when stitching them
            together. If zero the end time of each step overlaps the starting
            time of its following step. The default is 1e-3.

        """

        super().__init__()

        self._solns = soln
        self._sim = soln[0]._sim.copy()

        sv_size = self._sim._sv0.size

        self.message = []
        self.success = []
        self.status = []

        self.t = np.empty([0])
        self.y = np.empty([0, sv_size])
        self.yp = np.empty([0, sv_size])

        self.t_events = None
        self.y_events = None
        self.yp_events = None

        self.nfev = []
        self.njev = []

        self._timers = []

        for soln in self._solns:
            if self.t.size > 0:
                shift_t = self.t[-1] + soln.t + t_shift
            else:
                shift_t = soln.t

            if soln.t_events and self.t.size > 0:
                shift_t_events = self.t[-1] + soln.t_events + t_shift
            elif soln.t_events:
                shift_t_events = soln.t_events

            self.message.append(soln.message)
            self.success.append(soln.success)
            self.status.append(soln.status)

            self.t = np.hstack([self.t, shift_t])
            self.y = np.vstack([self.y, soln.y])
            self.yp = np.vstack([self.yp, soln.yp])

            if soln.t_events:
                if self.t_events is None:
                    self.t_events = shift_t_events
                    self.y_events = soln.y_events
                    self.yp_events = soln.yp_events
                else:
                    self.t_events = np.hstack([self.t_events, shift_t_events])
                    self.y_events = np.vstack([soln.y_events])
                    self.yp_events = np.vstack([soln.yp_events])

            self.nfev.append(soln.nfev)
            self.njev.append(soln.njev)

            self._timers.append(soln._timer)

        self._fill_vars()

    @property
    def solvetime(self) -> str:
        """
        Print a statement specifying how long IDASolver spent integrating.

        Returns
        -------
        solvetime : str
            An f-string with the total solver integration time in seconds.

        """
        return f"{sum(self._timers):.3f} s"

    def get_steps(self, idx: int | tuple) -> StepSolution | CycleSolution:
        """
        Return a subset of the solution.

        Parameters
        ----------
        idx : int | tuple
            The step index (int) or first/last indices (tuple) to return.

        Returns
        -------
        :class:`StepSolution` | :class:`CycleSolution`
            The returned solution subset. A StepSolution is returned if 'idx'
            is an int, and a CycleSolution will be returned for the range of
            requested steps when 'idx' is a tuple.

        """

        if isinstance(idx, int):
            return deepcopy(self._solns[idx])
        elif isinstance(idx, (tuple, list)):
            solns = self._solns[idx[0]:idx[1] + 1]
            return CycleSolution(*solns)
