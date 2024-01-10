:py:mod:`bmlite`
================

.. py:module:: bmlite

.. autoapi-nested-parse::

   .. rubric:: BatMods-lite

   Provides:

   1. Simple API to run pre-built battery models
   2. Library of kinetic and transport properties for common battery materials

   .. rubric:: How to use the documentation

   Documentation is accessible via python's ``help()`` function which prints
   docstrings from the specified package, module, function, class, etc. (e.g.,
   ``help(bmlite.SPM)``). In addition, you can access the documentation
   by calling the built-in ``bmlite.docs()`` method to open a locally run website.
   The website is recommended because it has search functionality and includes
   examples in jupyter notebooks.

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
   bmlite.format_ax



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
   array to integrate/solve the system. The ``inputs`` parameter is a ``tuple``
   that is used to pass any required user-defined ``*args`` to the function.

   :param residuals: Residual function like ``def residuals(t, y, yp, res, inputs)``.
                     For some examples, see :func:`bmlite.SPM.dae.residuals` and/or
                     :func:`bmlite.P2D.dae.residuals`.
   :type residuals: Callable
   :param \*\*kwargs: The keyword arguments specify the Sundials IDA solver options. A
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
   :type \*\*kwargs: dict, optional

   .. rubric:: Notes

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


.. py:function:: docs() -> None

   Opens a new tab in your browser with a locally run docs website.

   :returns: *None.*


.. py:function:: format_ax(ax: object) -> None

   Formats an ``axis`` object by adjusting the ticks.

   Specifically, the top and right ticks are added, minor ticks are turned on,
   and all ticks are set to face inward.

   :param ax: An ``axis`` instance from a ``matplotlib`` figure.
   :type ax: object

   :returns: *None.*


