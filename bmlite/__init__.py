"""
BatMods-lite
============

Provides:

1. Simple API to run pre-built battery models
2. Library of kinetic and transport properties for common battery materials

How to use the documentation
----------------------------
Documentation is accessible via python's ``help()`` function which prints
docstrings from the specified package, module, function, class, etc. (e.g.,
``help(bmlite.SPM)``). In addition, you can access the documentation
by calling the built-in ``bmlite.docs()`` method to open a locally run website.
The website is recommended because it has search functionality and includes
examples in jupyter notebooks.

Viewing documentation using IPython
-----------------------------------
Start IPython and import ``bmlite`` under the alias ``bm``. To see what's
available in ``bmlite``, type ``bm.<TAB>`` (where ``<TAB>`` refers to the TAB
key). To view the type hints and brief descriptions, type an open parenthesis
``(`` after any function, method, class, etc. (e.g., ``bm.Constants(``).
"""

from typing import Callable as _Callable

from numpy import ndarray as _ndarray
from scikits.odes import dae as _DAE

from . import P2D, SPM, materials


class Constants(object):
    __slots__ = []

    _F = 96485.3321e3
    _R = 8.3145e3

    def __init__(self) -> None:
        """
        Physical constants class with read-only attributes.
        """
        pass

    @property
    def F(self) -> float:
        """
        Faraday's constant [C/kmol].
        """
        return self._F

    @property
    def R(self) -> float:
        """
        Gas constant [J/kmol/K].
        """
        return self._R


class IDASolver(_DAE):
    """
    An IDA solver defined by a residuals function.

    For an ODE or DAE defined as ``M*y' = f(t, y)``, the residuals function is
    ``residuals = M*y' - f(t, y)``. This must be built as a python function
    with a signature like ``def resdiuals(t, y, yp, res, inputs) -> None``.
    The ``res`` parameter must be a 1D array the same size as ``y`` and ``yp``.
    Although the function returns ``None``, the solver uses the filled ``res``
    array to integrate/solve the system. The ``inputs`` parameter is a ``tuple``
    that is used to pass any required user-defined ``*args`` to the function.

    Parameters
    ----------
    residuals : Callable
        Residual function like ``def residuals(t, y, yp, res, inputs)``.
        For some examples, see :func:`bmlite.SPM.dae.residuals` and/or
        :func:`bmlite.P2D.dae.residuals`.

    **kwargs : dict, optional
        The keyword arguments specify the Sundials IDA solver options. A
        partial list of options/defaults is given below:

        ================== ==================================================
        Key                Description (type or options, default)
        ================== ==================================================
        rtol               relative tolerance (*float*, 1e-6)
        atol               absolute tolerance (*float*, 1e-12)
        user_data          the ``inputs`` parameter (*tuple*, ``None``)
        linsolver          linear solver (``{'dense', 'band'}``, ``'dense'``)
        lband              width of the lower band (*int*, 0)
        uband              width of the upper band (*int*, 0)
        max_step_size      max time step [s] (*float*, 0. -> unrestricted)
        rootfn             root/event function (*Callable*, ``None``)
        nr_rootfns         number of events in ``'rootfn'`` (*int*, ``0``)
        compute_initcond   vary y or yp to get a consistent initial condition
                           (``{'y0', 'yp0', None}``, ``'yp0'``)
        algebraic_vars_idx the indices of the algebraic variables in the
                           ``residuals`` y array (*list[int]*, ``None``)
        ================== ==================================================

    Notes
    -----
    * The solver name IDA stands for Implicit Differential-Algebraic solver. It
      is part of the `SUNDIALS`_ package, and is accessed here through the
      `scikits-odes`_ python wrapper.
    * The solver can be unstable if the ``algebraic_vars_idx`` keyword argument
      is not specified for DAEs.
    * The ``rootfn`` keyword argument must have the signature ``def f(t, y, yp,
      events, inputs) -> None`` where the ``events`` parameter is an array that
      is filled with root functions. If any element in ``events`` hits zero
      during the solver integration, the solver will exit.

      .. _SUNDIALS: https://sundials.readthedocs.io/
      .. _scikits-odes: https://bmcage.github.io/odes/dev/
    """

    def __init__(self, residuals, **kwargs) -> None:
        options = {}
        options["rtol"] = kwargs.get("rtol", 1e-6)
        options["atol"] = kwargs.get("atol", 1e-9)
        options["user_data"] = kwargs.get("user_data", None)
        options["linsolver"] = kwargs.get("linsolver", "dense")
        options["lband"] = kwargs.get("lband", 0)
        options["uband"] = kwargs.get("uband", 0)
        options["max_step_size"] = kwargs.get("max_step_size", 0.0)
        options["rootfn"] = kwargs.get("rootfn", None)
        options["nr_rootfns"] = kwargs.get("nr_rootfns", 0)

        options["old_api"] = kwargs.get("old_api", False)
        options["compute_initcond"] = kwargs.get("compute_initcond", "yp0")
        options["algebraic_vars_idx"] = kwargs.get("algebraic_vars_idx", None)

        super().__init__("ida", residuals, **options)


def docs() -> None:
    """
    Opens a new tab in your browser with a locally run docs website.

    Returns
    -------
    None.
    """

    import os
    import webbrowser

    sitepath = os.path.dirname(__file__) + "/docs/build/index.html"

    webbrowser.open_new_tab("file:" + sitepath)


def format_ax(ax: object) -> None:
    """
    Formats an ``axis`` object by adjusting the ticks.

    Specifically, the top and right ticks are added, minor ticks are turned on,
    and all ticks are set to face inward.

    Parameters
    ----------
    ax : object
        An ``axis`` instance from a ``matplotlib`` figure.

    Returns
    -------
    None.
    """

    from matplotlib.ticker import AutoMinorLocator

    if ax.get_xaxis().get_scale() != "log":
        ax.xaxis.set_minor_locator(AutoMinorLocator())

    if ax.get_yaxis().get_scale() != "log":
        ax.yaxis.set_minor_locator(AutoMinorLocator())

    ax.tick_params(axis="x", top=True, which="both", direction="in")
    ax.tick_params(axis="y", right=True, which="both", direction="in")


__version__ = "0.0.1"
