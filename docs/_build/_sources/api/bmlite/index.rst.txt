:py:mod:`bmlite`
================

.. py:module:: bmlite

.. autoapi-nested-parse::

   .. rubric:: BatMods-lite

   Battery Analysis and Training Models for Optimization and Degradation Studies
   (BATMODS) is a Python package with an API for pre-built battery models. The
   original purpose of the package was to quickly generate synthetic data for
   machine learning models to train with. However, the models are generally useful
   for any battery simulation or analysis. BATMODS-lite includes the following:

   1) A library and API for pre-built battery models
   2) Kinetic/transport properties for common battery materials

   .. rubric:: How to use the documentation

   Documentation is accessible via Python's ``help()`` function which prints
   docstrings from the specified package, module, function, class, etc. (e.g.,
   ``help(bmlite.SPM)``). In addition, you can access the documentation
   by calling the built-in ``bmlite.docs()`` method to open a locally run website.
   The website includes search functionality and examples, beyond the code
   docstrings.

   .. rubric:: Viewing documentation using IPython

   Start IPython and import ``bmlite`` under the alias ``bm``. To see what's
   available in ``bmlite``, type ``bm.<TAB>`` (where ``<TAB>`` refers to the TAB
   key). To view the type hints and brief descriptions, type an open parenthesis
   ``(`` after any function, method, class, etc. (e.g., ``bm.Constants(``).



Subpackages
-----------
.. toctree::
   :titlesonly:
   :maxdepth: 3

   P2D/index.rst
   SPM/index.rst
   materials/index.rst
   plotutils/index.rst


Submodules
----------
.. toctree::
   :titlesonly:
   :maxdepth: 1

   math/index.rst
   mesh/index.rst


Package Contents
----------------

Classes
~~~~~~~

.. autoapisummary::

   bmlite.Constants
   bmlite.IDASolver



Functions
~~~~~~~~~

.. autoapisummary::

   bmlite.docs



.. py:class:: Constants




   
   Physical constants class with read-only attributes.

   .. py:property:: F
      :type: float

      Faraday's constant [C/kmol].

   .. py:property:: R
      :type: float

      Gas constant [J/kmol/K].


.. py:class:: IDASolver(residuals, **kwargs)




   An IDA solver defined by a residuals function.

   For an ODE or DAE defined as ``M*y' = f(t, y)``, the residuals function is
   ``residuals = M*y' - f(t, y)``. This must be built as a python function
   with a signature like ``def resdiuals(t, y, yp, res, inputs) -> None``.
   The ``res`` parameter must be a 1D array the same size as ``y`` and ``yp``.
   Although the function returns ``None``, the solver uses the filled ``res``
   array to integrate/solve the system. The ``inputs`` parameter is a *tuple*
   that is used to pass any required user-defined ``*args`` to the function.

   :param residuals: Residuals function ``def residuals(t, y, yp, res, inputs) -> None``.
                     For some examples, see :func:`bmlite.SPM.dae.residuals` and/or
                     :func:`bmlite.P2D.dae.residuals`.
   :type residuals: Callable
   :param \*\*kwargs: The keyword arguments specify the Sundials IDA solver options. A
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
                      algidx       algebraic variable indices in y (*list[int]*, ``None``)
                      max_t_step   maximum time step [s] (*float*, 0. -> unrestricted)
                      ============ =========================================================
   :type \*\*kwargs: dict, optional

   .. rubric:: Notes

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

   .. py:method:: init_step(t0, y0, yp0)

      Solve the ODE/DAE for a consistent initial condition at ``t = t0``.

      :param t0: Initial time [s].
      :type t0: float
      :param y0: Array of state variables at ``t = t0``.
      :type y0: 1D array
      :param yp0: Array of state variable time derivatives at ``t = t0``.
      :type yp0: 1D array

      :returns: **idasol** (*NamedTuple*) -- Solution returned by SUNDIALS IDA integrator at ``t = t0``.


   .. py:method:: solve(t_span, y0, yp0)

      Solve the ODE/DAE system and save the solution at each time in
      ``t_span``.

      :param t_span: Array of times [s] to store the solution. ``t_span[0]`` must be
                     the start time and ``t_span[-1]`` the final time.
      :type t_span: 1D array
      :param y0: Array of state variables at ``t = t_span[0]``.
      :type y0: 1D array
      :param yp0: Array of state variable time derivatives at ``t = t_span[0]``.
      :type yp0: 1D array

      :returns: **idasol** (*NamedTuple*) -- Solution returned by SUNDIALS IDA integrator over ``t_span``.



.. py:function:: docs()

   Opens a new tab in your browser with a locally run docs website.

   :returns: *None.*


