bmlite
======

.. py:module:: bmlite

.. autoapi-nested-parse::

   .. rubric:: BATMODS-lite

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
   :maxdepth: 1

   /api/bmlite/P2D/index
   /api/bmlite/SPM/index
   /api/bmlite/materials/index
   /api/bmlite/plotutils/index


Submodules
----------

.. toctree::
   :maxdepth: 1

   /api/bmlite/mathutils/index
   /api/bmlite/mesh/index


Classes
-------

.. autoapisummary::

   bmlite.Constants
   bmlite.Experiment
   bmlite.IDAResult
   bmlite.IDASolver


Functions
---------

.. autoapisummary::

   bmlite.templates


Package Contents
----------------

.. py:class:: Constants

   Physical constants.


   .. py:property:: F
      :type: float

      Faraday's constant [C/kmol].


   .. py:property:: R
      :type: float

      Gas constant [J/kmol/K].


.. py:class:: Experiment(**kwargs)

   Experiment builder.

   A class to define an experimental protocol. Use the add_step() method
   to add a series of sequential steps. Each step defines a control mode,
   a constant or time-dependent load profile, a time span, and optional
   limiting criteria to stop the step early if a specified event/state is
   detected.

   :param \*\*kwargs: IDASolver keyword arguments that span all steps.
   :type \*\*kwargs: dict, optional

   .. seealso::

      :obj:`~bmlite.IDASolver`
          The solver class, with documentation for most keyword arguments that you might want to adjust.


   .. py:method:: add_step(mode, value, tspan, limits = None, **kwargs)

      Add a step to the experiment.

      :param mode: Control mode, {'current_A', 'current_C', 'voltage_V', 'power_W'}.
      :type mode: str
      :param value: Value of boundary contion mode, in the appropriate units.
      :type value: float | Callable
      :param tspan: Relative times for recording solution [s]. Providing a tuple as
                    (t_max: float, Nt: int) or (t_max: float, dt: float) constructs
                    tspan using ``np.linspace`` or ``np.arange``, respectively. Given
                    an array uses the values supplied as the evaluation times. Arrays
                    must be monotonically increasing and start with zero. See the notes
                    for more information.
      :type tspan: tuple | 1D np.array
      :param limits: Stopping criteria for the new step, must be entered in sequential
                     name/value pairs. Allowable names are {'current_A', 'current_C',
                     'voltage_V', 'power_W', 'time_s', 'time_min', 'time_h'}. Values for
                     each limit should immediately follow a corresponding name and match
                     its units. Time limits are in reference to total experiment time.
                     The default is None.
      :type limits: tuple[str, float], optional
      :param \*\*kwargs: IDASolver keyword arguments specific to the new step only.
      :type \*\*kwargs: dict, optional

      :returns: *None.*

      :raises ValueError: 'mode' is invalid.
      :raises ValueError: A 'limits' name is invalid.
      :raises ValueError: 'tspan' tuple must be length 2.
      :raises TypeError: 'tspan[1]' must be type int or float.
      :raises ValueError: 'tspan' arrays must be one-dimensional.
      :raises ValueError: 'tspan[0]' must be zero when given an array.
      :raises ValueError: 'tspan' array length must be at least two.
      :raises ValueError: 'tspan' arrays must be monotonically increasing.

      .. seealso::

         :obj:`~bmlite.IDASolver`
             The solver class, with documentation for most keyword arguments that you might want to adjust.

      .. rubric:: Notes

      For time-dependent loads, use a Callable for 'value' with a function
      signature like ``def load(t: float) -> float``, where 't' is the step's
      relative time, in seconds.

      Solution times are constructed and saved depending on the 'tspan' input
      types that were supplied:

      * Given (float, int):
          ``tspan = np.linspace(0., tspan[0], tspan[1])``
      * Given (float, float):
          ``tspan = np.arange(0., tspan[0], tspan[1])``

          In this case, 't_max' is also appended to the end. This results
          in the final 'dt' being different from the others if 't_max' is
          not evenly divisible by the given 'dt'.
      * Given 1D np.array:
          When you provide a numpy array it is checked for compatibility.
          If the array is not 1D, is not monotonically increasing, or starts
          with a value other than zero then an error is raised.



   .. py:method:: print_steps()

      Prints a formatted/readable list of steps.

      :returns: *None.*



   .. py:property:: num_steps
      :type: int

      Return number of steps.

      :returns: **num_steps** (*int*) -- Number of steps.


   .. py:property:: steps
      :type: list[dict]

      Return steps list.

      :returns: **steps** (*list[dict]*) -- List of the step dictionaries.


.. py:class:: IDAResult(**kwargs)



   Results class for IDA solver.

   Inherits from :class:`~sksundae.common.RichResult`. The solution class
   groups output from :class:`IDA` into an object with the fields:

   :param message: Human-readable description of the status value.
   :type message: str
   :param success: True if the solver was successful (status >= 0). False otherwise.
   :type success: bool
   :param status: Reason for the algorithm termination. Negative values correspond
                  to errors, and non-negative values to different successful criteria.
   :type status: int
   :param t: Solution time(s). The dimension depends on the method. Stepwise
             solutions will only have 1 value whereas solutions across a full
             'tspan' will have many.
   :type t: ndarray, shape(n,)
   :param y: State variable values at each solution time. Rows correspond to
             indices in 't' and columns match indexing from 'y0'.
   :type y: ndarray, shape(n, m)
   :param yp: State variable time derivate values at each solution time. Row
              and column indexing matches 'y'.
   :type yp: ndarray, shape(n, m)
   :param i_events: Provides an array for each detected event 'k' specifying indices
                    for which event(s) occurred. ``i_events[k,i] != 0`` if 'events[i]'
                    occurred at 't_events[k]'. The sign of 'i_events' indicates the
                    direction of zero-crossing:

                        * -1 indicates 'events[i]' was decreasing
                        * +1 indicates 'events[i]' was increasing

                    Output for 'i_events' will be None when either 'eventsfn' was None
                    or if no events occurred during the solve.
   :type i_events: ndarray, shape(k, num_events) or None
   :param t_events: Times at which events occurred or None if 'eventsfn' was None or
                    no events were triggered during the solve.
   :type t_events: ndarray, shape(k,) or None
   :param y_events: State variable values at each 't_events' value or None. Rows and
                    columns correspond to 't_events' and 'y0' indexing, respectively.
   :type y_events: ndarray, shape(k, m) or None
   :param yp_events: State variable time derivative values at each 't_events' value or
                     None. Row and column indexing matches 'y_events'.
   :type yp_events: ndarray, shape(k, m) or None
   :param nfev: Number of times that 'resfn' was evaluated.
   :type nfev: int
   :param njev: Number of times the Jacobian was evaluated, 'jacfn' or internal
                finite difference method.
   :type njev: int

   .. rubric:: Notes

   Terminal events are appended to the end of 't', 'y', and 'yp'. However,
   if an event was not terminal then it will only appear in '\*_events'
   outputs and not within the main output arrays.

   'nfev' and 'njev' are cumulative for stepwise solution approaches. The
   values are reset each time 'init_step' is called.


.. py:class:: IDASolver(resfn, **options)



   SUNDIALS IDA solver.

   This class wraps the implicit differential algebraic (IDA) solver from
   SUNDIALS. It can be used to solve both ordinary differential equations
   (ODEs) and differiential agebraic equatinos (DAEs).

   :param resfn: Residual function with signature ``f(t, y, yp, res[, userdata])``.
                 If 'resfn' has return values, they are ignored. Instead of using
                 returns, the solver interacts directly with the 'res' array memory.
                 For more info see the notes.
   :type resfn: Callable
   :param \*\*options: Keyword arguments to describe the solver options. A full list of
                       names, types, descriptions, and defaults is given below.
   :type \*\*options: dict, optional
   :param userdata: Additional data object to supply to all user-defined callables. If
                    'resfn' takes in 5 arguments, including the optional 'userdata',
                    then this option cannot be None (default). See notes for more info.
   :type userdata: object or None, optional
   :param calc_initcond: Specifies which initial condition, if any, to calculate prior to
                         the first time step. The options 'y0' and 'yp0' will correct 'y0'
                         or 'yp0' values at 't0', respectively. When not None (default),
                         the 'calc_init_dt' value should be used to specify the direction
                         of integration.
   :type calc_initcond: {'y0', 'yp0', None}, optional
   :param calc_init_dt: Relative time step to take during the initial condition correction.
                        Positive vs. negative values provide the direction of integration
                        as forwards or backwards, respectively. The default is 0.01.
   :type calc_init_dt: float, optional
   :param algebraic_idx: Specifies indices 'i' in the 'y[i]' state variable array that are
                         purely algebraic. This option should always be provided for DAEs;
                         otherwise, the solver can be unstable. The default is None.
   :type algebraic_idx: array_like[int] or None, optional
   :param first_step: Specifies the initial step size. The default is 0, which uses an
                      estimated value internally determined by SUNDIALS.
   :type first_step: float, optional
   :param min_step: Minimum allowable step size. The default is 0.
   :type min_step: float, optional
   :param max_step: Maximum allowable step size. Use 0 (default) for unbounded steps.
   :type max_step: float, optional
   :param rtol: Relative tolerance. For example, 1e-4 means errors are controlled
                to within 0.01%. It is recommended to not use values larger than
                1e-3 nor smaller than 1e-15. The default is 1e-5.
   :type rtol: float, optional
   :param atol: Absolute tolerance. Can be a scalar float to apply the same value
                for all state variables, or an array with a length matching 'y' to
                provide tolerances specific to each variable. The default is 1e-6.
   :type atol: float or array_like[float], optional
   :param linsolver: Choice of linear solver. When using 'band', don't forget to provide
                     'lband' and 'uband' values. The default is 'dense'.
   :type linsolver: {'dense', 'band'}, optional
   :param lband: Lower Jacobian bandwidth. Given a DAE system ``0 = F(t, y, yp)``,
                 the Jacobian is ``J = dF_i/dy_j + c_j*dF_i/dyp_j`` where 'c_j' is
                 determined internally based on both step size and order. 'lband'
                 should be set to the max distance between the main diagonal and the
                 non-zero elements below the diagonal. This option cannot be None
                 (default) if 'linsolver' is 'band'. Use zero if no values are below
                 the main diagonal.
   :type lband: int or None, optional
   :param uband: Upper Jacobian bandwidth. See 'lband' for the Jacobian description.
                 'uband' should be set to the max distance between the main diagonal
                 and the non-zero elements above the diagonal. This option cannot be
                 None (default) if 'linsolver' is 'band'. Use zero if no elements
                 are above the main diagonal.
   :type uband: int or None, optional
   :param max_order: Specifies the maximum order for the linear multistep BDF method.
                     The value must be in the range [1, 5]. The default is 5.
   :type max_order: int, optional
   :param max_num_steps: Specifies the maximum number of steps taken by the solver in each
                         attempt to reach the next output time. The default is 500.
   :type max_num_steps: int, optional
   :param max_nonlin_iters: Specifies the maximum number of nonlinear solver iterations in one
                            step. The default is 4.
   :type max_nonlin_iters: int, optional
   :param max_conv_fails: Specifies the max number of nonlinear solver convergence failures
                          in one step. The default is 10.
   :type max_conv_fails: int, optional
   :param constraints_idx: Specifies indices 'i' in the 'y' state variable array for which
                           inequality constraints should be applied. Constraints types must be
                           specified in 'constraints_type', see below. The default is None.
   :type constraints_idx: array_like[int] or None, optional
   :param constraints_type: If 'constraints_idx' is not None, then this option must include an
                            array of equal length specifying the types of constraints to apply.
                            Values should be in ``{-2, -1, 1, 2}`` which apply ``y[i] < 0``,
                            ``y[i] <= 0``, ``y[i] >=0,`` and ``y[i] > 0``, respectively. The
                            default is None.
   :type constraints_type: array_like[int] or None, optional
   :param eventsfn: Events function with signature ``g(t, y, yp, events[, userdata])``.
                    Return values from this function are ignored. Instead, the solver
                    directly interacts with the 'events' array. Each 'events[i]' should
                    be an expression that triggers an event when equal to zero. If None
                    (default), no events are tracked. See the notes for more info.

                    The 'num_events' option is required when 'eventsfn' is not None so
                    memory can be allocated for the events array. The events function
                    can also have the following attributes:

                        terminal: list[bool, int], optional
                            A list with length 'num_events' that tells how the solver
                            how to respond to each event. If boolean, the solver will
                            terminate when True and will simply record the event when
                            False. If integer, termination occurs at the given number
                            of occurrences. The default is ``[True]*num_events``.
                        direction: list[int], optional
                            A list with length 'num_events' that tells the solver which
                            event directions to track. Values must be in ``{-1, 0, 1}``.
                            Negative values will only trigger events when the slope is
                            negative (i.e., 'events[i]' went from positive to negative).
                            Alternatively, positive values track events with positive
                            slope. If zero, either direction triggers the event. When
                            not assigned, ``direction = [0]*num_events``.

                    You can assign attributes like ``eventsfn.terminal = [True]`` to
                    any function in Python, after it has been defined.
   :type eventsfn: Callable or None, optional
   :param num_events: Number of events to track. Must be greater than zero if 'eventsfn'
                      is not None. The default is 0.
   :type num_events: int, optional
   :param jacfn: Jacobian function like ``J(t, y, yp, res, cj, JJ[, userdata])``.
                 The function should fill the pre-allocated 2D matrix 'JJ' with the
                 values defined by ``JJ[i,j] = dres_i/dy_j + cj*dres_i/dyp_j``. An
                 internal finite difference method is applied when None (default).
                 As with other user-defined callables, return values from 'jacfn'
                 are ignored. See notes for more info.
   :type jacfn: Callable or None, optional

   .. rubric:: Notes

   Return values from 'resfn', 'eventsfn', and 'jacfn' are ignored by the
   solver. Instead the solver directly reads from pre-allocated memory.
   The 'res', 'events', and 'JJ' arrays from each user-defined callable
   should be filled within each respective function. When setting values
   across the entire array/matrix at once, don't forget to use ``[:]`` to
   fill the existing array rather than overwriting it. For example, using
   ``res[:] = F(t, y, yp)`` is correct whereas ``res = F(t, y, yp)`` is
   not. Using this method of pre-allocated memory helps pass data between
   Python and the SUNDIALS C functions. It also keeps the solver fast,
   especially for large problems.

   When 'resfn' (or 'eventsfn', or 'jacfn') require data outside of their
   normal arguments, you can supply 'userdata' as an option. When given,
   'userdata' must appear in the function signatures for ALL of 'resfn',
   'eventsfn' (when not None), and 'jacfn' (when not None), even if it is
   not used in all of these functions. Note that 'userdata' only takes up
   one argument position; however, 'userdata' can be any Python object.
   Therefore, to pass more than one extra argument you should pack all of
   the data into a single tuple, dict, dataclass, etc. and pass them all
   together as 'userdata'. The data can be unpacked as needed within a
   function.

   .. rubric:: Examples

   The following example solves the Robertson problem, which is a classic
   test problem for programs that solve stiff ODEs. A full description of
   the problem is provided by `MATLAB`_. Note that while initializing the
   solver, ``algebraic_idx=[2]`` specifies ``y[2]`` is purely algebraic,
   and ``calc_initcond='yp0'`` tells the solver to determine the values
   for 'yp0' at 'tspan[0]' before starting to integrate. That is why 'yp0'
   can be initialized as an array of zeros even though plugging in 'y0'
   to the residuals expressions actually gives ``yp0 = [-0.04, 0.04, 0]``.
   The initialization is checked against the correct answer after solving.

   .. _MATLAB:
       https://mathworks.com/help/matlab/math/
       solve-differential-algebraic-equations-daes.html

   .. code-block:: python

       import numpy as np
       import sksundae as sun
       import matplotlib.pyplot as plt

       def resfn(t, y, yp, res):
           res[0] = yp[0] + 0.04*y[0] - 1e4*y[1]*y[2]
           res[1] = yp[1] - 0.04*y[0] + 1e4*y[1]*y[2] + 3e7*y[1]**2
           res[2] = y[0] + y[1] + y[2] - 1.0

       solver = sun.ida.IDA(resfn, algebraic_idx=[2], calc_initcond='yp0')

       tspan = np.hstack([0, 4*np.logspace(-6, 6)])
       y0 = np.array([1, 0, 0])
       yp0 = np.zeros_like(y0)

       soln = solver.solve(tspan, y0, yp0)
       assert np.allclose(soln.yp[0], [-0.04, 0.04, 0], rtol=1e-3)

       soln.y[:, 1] *= 1e4  # scale y[1] so it is visible in the figure
       plt.semilogx(soln.t, soln.y)
       plt.show()


.. py:function:: templates(model, file = None)

   Print simulation templates.

   If no input, a list of available templates will print. Otherwise, given
   a file name or index, a template will print.

   :param model: Model package name, e.g., 'SPM'. Input is case insensitive. All model
                 names are converted to uppercase within this function.
   :type model: str
   :param file: File name or index. The default is None.
   :type file: str | int, optional

   :returns: *None.*

   :raises AttributeError: 'model' is not a valid model package.
   :raises FileNotFoundError: 'model' has no 'templates' directory.


