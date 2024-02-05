"""
BatMods-lite
============

Provides:

1. A library and API for pre-built battery models
2. Kinetic/transport properties for common battery materials

How to use the documentation
----------------------------
Documentation is accessible via Python's ``help()`` function which prints
docstrings from the specified package, module, function, class, etc. (e.g.,
``help(bmlite.SPM)``). In addition, you can access the documentation
by calling the built-in ``bmlite.docs()`` method to open a locally run website.
The website includes search functionality and examples, beyond the code
docstrings.

Viewing documentation using IPython
-----------------------------------
Start IPython and import ``bmlite`` under the alias ``bm``. To see what's
available in ``bmlite``, type ``bm.<TAB>`` (where ``<TAB>`` refers to the TAB
key). To view the type hints and brief descriptions, type an open parenthesis
``(`` after any function, method, class, etc. (e.g., ``bm.Constants(``).
"""


from bmlite.version import __version__
from numpy import ndarray as _ndarray
from scikits.odes import dae as _DAE

from . import materials
from . import P2D
from . import plotutils
from . import SPM

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


class IDASolver(object):
    """
    An IDA solver defined by a residuals function.

    For an ODE or DAE defined as ``M*y' = f(t, y)``, the residuals function is
    ``residuals = M*y' - f(t, y)``. This must be built as a python function
    with a signature like ``def resdiuals(t, y, yp, res, inputs) -> None``.
    The ``res`` parameter must be a 1D array the same size as ``y`` and ``yp``.
    Although the function returns ``None``, the solver uses the filled ``res``
    array to integrate/solve the system. The ``inputs`` parameter is a *tuple*
    that is used to pass any required user-defined ``*args`` to the function.

    Parameters
    ----------
    residuals : Callable
        Residuals function ``def residuals(t, y, yp, res, inputs) -> None``.
        For some examples, see :func:`bmlite.SPM.dae.residuals` and/or
        :func:`bmlite.P2D.dae.residuals`.

    **kwargs : dict, optional
        The keyword arguments specify the Sundials IDA solver options. A
        partial list of options/defaults is given below:

        ============ =========================================================
        Key          Description (*type* or {options}, default)
        ============ =========================================================
        rtol         relative tolerance (*float*, 1e-6)
        atol         absolute tolerance (*float*, 1e-12)
        user_data    the ``inputs`` parameter (*tuple*, ``None``)
        linsolver    linear solver (``{'dense', 'band'}``, ``'dense'``)
        lband        width of the lower band (*int*, 0)
        uband        width of the upper band (*int*, 0)
        rootfn       root/event function (*Callable*, ``None``)
        nr_rootfns   number of events in ``'rootfn'`` (*int*, ``0``)
        initcond     unknown variable set (``{'y0', 'yp0', None}``, ``'yp0'``)
        algidx       algebraic variable indicies in y (*list[int]*, ``None``)
        max_t_step   maximum time step [s] (*float*, 0. -> unrestricted)
        ============ =========================================================

    Notes
    -----
    * The solver name IDA stands for Implicit Differential-Algebraic solver.
      It is part of the `SUNDIALS`_ package, and is accessed here through the
      `scikits-odes`_ python wrapper.
    * The solver can be unstable if the ``algidx`` keyword argument is not
      specified for DAEs.
    * The ``rootfn`` keyword argument must have the signature ``def f(t, y, yp,
      events, inputs) -> None`` where the ``events`` parameter is an array that
      is filled with root functions. If any element in ``events`` hits zero
      during the solver integration, the solver will exit.

      .. _SUNDIALS: https://sundials.readthedocs.io/
      .. _scikits-odes: https://bmcage.github.io/odes/dev/
    """

    __slots__ = ['_integrator']

    def __init__(self, residuals, **kwargs) -> None:

        # Overwrite scikits.odes defaults w/ some keys renamed
        options = {}
        options['rtol'] = kwargs.pop('rtol', 1e-6)
        options['atol'] = kwargs.pop('atol', 1e-9)
        options['user_data'] = kwargs.pop('user_data', None)

        options['linsolver'] = kwargs.pop('linsolver', 'dense')
        options['lband'] = kwargs.pop('lband', 0)
        options['uband'] = kwargs.pop('uband', 0)

        options['rootfn'] = kwargs.pop('rootfn', None)
        options['nr_rootfns'] = kwargs.pop('nr_rootfns', 0)

        options['compute_initcond'] = kwargs.pop('initcond', 'yp0')
        options['algebraic_vars_idx'] = kwargs.pop('algidx', None)
        options['max_step_size'] = kwargs.pop('max_t_step', 0.)

        options['old_api'] = False

        # Collect new defaults and any extra user kwargs
        options = {**options, **kwargs}

        _DAE.__init__(self, 'ida', residuals, **options)

    def solve(self, t_span: _ndarray, y0: _ndarray, yp0: _ndarray) -> object:
        """
        Solve the ODE/DAE system and save the solution at each time in
        ``t_span``.

        Parameters
        ----------
        t_span : 1D array
            Array of times [s] to store the solution. ``t_span[0]`` must be
            the start time and ``t_span[-1]`` the final time.

        y0 : 1D array
            Array of state variables at ``t = t_span[0]``.

        yp0 : 1D array
            Array of state variable time derivatives at ``t = t_span[0]``.

        Returns
        -------
        idasol : NamedTuple
            Solution returned by SUNDIALS IDA integrator over ``t_span``.
        """

        idasol = _DAE.solve(self, t_span, y0, yp0)

        return idasol

    def init_step(self, t0: float, y0: _ndarray, yp0: _ndarray) -> object:
        """
        Solve the ODE/DAE for a consistent initial condition at ``t = t0``.

        Parameters
        ----------
        t0 : float
            Initial time [s].

        y0 : 1D array
            Array of state variables at ``t = t0``.

        yp0 : 1D array
            Array of state variable time derivatives at ``t = t0``.

        Returns
        -------
        idasol : NamedTuple
            Solution returned by SUNDIALS IDA integrator at ``t = t0``.
        """

        idasol = _DAE.init_step(self, t0, y0, yp0)

        return idasol


def docs() -> None:
    """
    Opens a new tab in your browser with a locally run docs website.

    Returns
    -------
    None.
    """

    import os, webbrowser

    sitepath = os.path.dirname(__file__) + '/../docs/_build/index.html'

    webbrowser.open_new_tab('file://' + os.path.realpath(sitepath))


def _templates(model__file__: str, model_name: str, sim: str | int = None,
               exp: str | int = None) -> None:
    """
    Print simulation and/or experiment templates. If both ``sim`` and ``exp``
    are ``None``, a list of available templates will be printed. Otherwise,
    if a name or index is given, that template will print to the console.

    Parameters
    ----------
    model__file__ : str
        The module ``__file__`` attribute. This sets the local path to make
        sure templates are pulled from the correct model path.

    model_name : str
        Name for the model package. This ensures template lists are labeled
        correctly.

    sim : str | int, optional
        Simulation template file name or index. The default is ``None``.

    exp : str | int, optional
        Experiment template file name or index. The default is ``None``.

    Returns
    -------
    None.
    """

    import os, json
    from pathlib import Path

    from ruamel.yaml import YAML

    dirname = os.path.dirname(model__file__)

    simlist = os.listdir(dirname + '/default_sims/')
    explist = os.listdir(dirname + '/default_exps/')

    if sim is None and exp is None:
        print('\nTemplates for ' + model_name + ' simulations:')
        for i, file in enumerate(simlist):
            print('- [' + str(i) + '] ' + file.removesuffix('.yaml'))

        print('\nTemplates for ' + model_name + ' experiments:')
        for i, file in enumerate(explist):
            print('- [' + str(i) + '] ' + file.removesuffix('.yaml'))

    if isinstance(sim, str):
        if '.yaml' not in sim:
            sim += '.yaml'

        simfile = sim
    elif isinstance(sim, int):
        simfile = simlist[sim]

    if sim is not None:
        print('\n' + '=' * 30 + '\n' + simfile + '\n' + '=' * 30)
        with open(dirname + '/default_sims/' + simfile, 'r') as f:
            print('\n' + f.read())

    yaml = YAML()

    if isinstance(exp, str):
        if '.yaml' not in exp:
            exp += '.yaml'

        expfile = exp
    elif isinstance(exp, int):
        expfile = explist[exp]

    if exp is not None:
        print('\n' + '=' * 30 + '\n' + expfile + '\n' + '=' * 30)
        expdict = yaml.load(Path(dirname + '/default_exps/' + expfile))
        print('exp = ' + json.dumps(expdict, indent=4))


__version__ = '0.0.1'
