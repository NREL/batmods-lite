:py:mod:`bmlite.SPM.solutions`
==============================

.. py:module:: bmlite.SPM.solutions

.. autoapi-nested-parse::

   .. rubric:: Solution Classes

   This package contains class definitions for each experimental solution. Most
   solutions inherit from :class:`~bmlite.SPM.solutions.BaseSolution`, however,
   their verification processes are unique.



Package Contents
----------------

Classes
~~~~~~~

.. autoapisummary::

   bmlite.SPM.solutions.BaseSolution
   bmlite.SPM.solutions.CCSolution
   bmlite.SPM.solutions.CPSolution
   bmlite.SPM.solutions.CVSolution




.. py:class:: BaseSolution(sim, exp)




   Base methods for all P2D Solution classes.

   When initialized, a copy of the SPM simulation object and experimental
   details for the solution instance are stored, all other instance
   attributes are set to ``None``. These are filled later with a "fill"
   method (e.g., ``ida_fill()`` or ``dict_fill()``).

   :param sim: The SPM Simulation instance used to produce the solution.
   :type sim: SPM Simulation object
   :param exp: Experiment dictionary. Specific key/value pairs are dependent on
               the experiment that was run.
   :type exp: dict

   .. py:property:: classname
      :type: str

      Class name.

      :returns: **classname** (*str*) -- Name of current class.

   .. py:property:: message
      :type: str

      Solver exit message.

      :returns: **message** (*str*) -- Exit message from Sundials IDA solver.

   .. py:property:: onroot
      :type: bool

      Solver root(event) exit status.

      :returns: **onroot** (*bool*) -- ``True`` if exit on root/event, ``False`` otherwise.

   .. py:property:: success
      :type: bool

      Solver overall exit status.

      :returns: **success** (*bool*) -- ``True`` if no errors, ``False`` otherwise.

   .. py:property:: t
      :type: numpy.ndarray

      Saved solution times [s].

      :returns: **t** (*1D array*) -- Solution times [s] returned by solver.

   .. py:property:: y
      :type: numpy.ndarray

      Saved spatiotemporal solution variables [units].

      :returns: **y** (*2D array*) -- Solution variables [units] returned by solver.

   .. py:property:: ydot
      :type: numpy.ndarray

      Saved solution variable time derivates (i.e., dy/dt) [units].

      :returns: **ydot** (*2D array*) -- Solution variable time derivatives [units] returned by solver.

   .. py:method:: dict_fill(sol)

      Fill the instance attributes using a dictionary.

      :param sol: Solution dictionary with the following key/value pairs:

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
      :type sol: dict

      :returns: *None.*


   .. py:method:: ida_fill(sol, solvetime)

      Fill the instance attributes using the SolverReturn object from the
      Sundials IDA solver.

      :param sol: Sundials IDA SolverReturn object with the following attributes:

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
      :type sol: object
      :param solvetime: Solver integration time, in seconds.
      :type solvetime: float

      :returns: *None.*


   .. py:method:: plot(*args)

      Generates requested plots based on ``*args``.

      :param \*args: Use any number of the following arguments to see the described
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
      :type \*args: str

      :returns: *None.*


   .. py:method:: post()


   .. py:method:: report()

      Prints the experiment details and solution success report.

      :returns: *None.*


   .. py:method:: save_sliced(savename, overwrite = False)

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

      :param savename: Either a file name or the absolute/relative file path. The ``.npz``
                       extension will be added to the end of the string if it is not
                       already there. If only the file name is given, the file will be
                       saved in the user's current working directory.
      :type savename: str
      :param overwrite: A flag to overwrite an existing ``.npz`` file with the same name
                        if one exists. The default is ``False``.
      :type overwrite: bool, optional

      :returns: *None.*


   .. py:method:: solvetime(units = 's')

      Print solve time (not including pre/post processing).

      :param units: Units for printed time ``{'s', 'min', 'h'}``. The default is
                    ``'s'``.
      :type units: str, optional

      :returns: **solvetime** (*str*) -- Time for Sundials IDA to solve problem in [units].


   .. py:method:: to_dict()

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

      :param None.:

      :returns: **sol_dict** (*dict*) -- A dictionary containing the solution



.. py:class:: CCSolution(sim, exp)




   Constant current solution for SPM simulations.

   Base: :class:`~bmlite.SPM.solutions.BaseSolution`

   When initialized, a copy of the SPM simulation object and experimental
   details for the solution instance are stored, all other instance
   attributes are set to ``None``. These are filled later with a "fill"
   method (e.g., ``ida_fill()`` or ``dict_fill()``).

   :param sim: The SPM Simulation instance used to produce the solution.
   :type sim: SPM Simulation object
   :param exp: Experiment dictionary. Specific key/value pairs are dependent on
               the experiment that was run.
   :type exp: dict

   .. py:property:: classname
      :type: str

      Class name. Overwrites ``BaseSolution``.

      :returns: **classname** (*str*) -- Name of current class.

   .. py:method:: verify(plotflag = False, rtol = 0.005, atol = 0.001)

      Verifies the solution is mathematically consistent.

      Specifically, for a constant current experiment, this method checks
      that the calculated current was within tolerance of the boundary
      condition. In addition, this method checks the anodic and cathodic
      Faradaic currents against the current at each time step and the
      solid-phase lithium conservation.

      If the verification returns ``False``, you can see which checks failed
      using ``plotflag``. Any subplots shaded grey failed their test.
      Failures generally suggest that the solver's relative and/or absolute
      tolerance should be adjusted, and/or that the discretization should be
      modified to adjust the mesh.

      :param plotflag: A flag to see plots showing the verification calculations. The
                       default is ``False``.
      :type plotflag: bool, optional
      :param rtol: Relative tolerance for comparisons. The default is 5e-3.
      :type rtol: float, optional
      :param atol: Absolute tolerance for comparisons. The default is 1e-3.
      :type atol: float, optional

      :returns: **checks** (*bool*) -- ``True`` is all checks are satisfied, ``False`` otherwise.



.. py:class:: CPSolution(sim, exp)




   Constant power solution for SPM simulations.

   Base: :class:`~bmlite.SPM.solutions.BaseSolution`

   When initialized, a copy of the SPM simulation object and experimental
   details for the solution instance are stored, all other instance
   attributes are set to ``None``. These are filled later with a "fill"
   method (e.g., ``ida_fill()`` or ``dict_fill()``).

   :param sim: The SPM Simulation instance used to produce the solution.
   :type sim: SPM Simulation object
   :param exp: Experiment dictionary. Specific key/value pairs are dependent on
               the experiment that was run.
   :type exp: dict

   .. py:property:: classname
      :type: str

      Class name. Overwrites ``BaseSolution``.

      :returns: **classname** (*str*) -- Name of current class.

   .. py:method:: verify(plotflag = False, rtol = 0.005, atol = 0.001)

      Verifies the solution is mathematically consistent.

      Specifically, for a constant power experiment, this method checks
      that the calculated power was within tolerance of the boundary
      condition. In addition, this method checks the anodic and cathodic
      Faradaic currents against the current at each time step and the
      solid-phase lithium conservation.

      If the verification returns ``False``, you can see which checks failed
      using ``plotflag``. Any subplots shaded grey failed their test.
      Failures generally suggest that the solver's relative and/or absolute
      tolerance should be adjusted, and/or that the discretization should be
      modified to adjust the mesh.

      :param plotflag: A flag to see plots showing the verification calculations. The
                       default is ``False``.
      :type plotflag: bool, optional
      :param rtol: Relative tolerance for comparisons. The default is 5e-3.
      :type rtol: float, optional
      :param atol: Absolute tolerance for comparisons. The default is 1e-3.
      :type atol: float, optional

      :returns: **checks** (*bool*) -- ``True`` is all checks are satisfied, ``False`` otherwise.



.. py:class:: CVSolution(sim, exp)




   Constant voltage solution for SPM simulations.

   Base: :class:`~bmlite.SPM.solutions.BaseSolution`

   When initialized, a copy of the SPM simulation object and experimental
   details for the solution instance are stored, all other instance
   attributes are set to ``None``. These are filled later with a "fill"
   method (e.g., ``ida_fill()`` or ``dict_fill()``).

   :param sim: The SPM Simulation instance used to produce the solution.
   :type sim: SPM Simulation object
   :param exp: Experiment dictionary. Specific key/value pairs are dependent on
               the experiment that was run.
   :type exp: dict

   .. py:property:: classname
      :type: str

      Class name. Overwrites ``BaseSolution``.

      :returns: **classname** (*str*) -- Name of current class.

   .. py:method:: verify(plotflag = False, rtol = 0.005, atol = 0.001)

      Verifies the solution is mathematically consistent.

      Specifically, for a constant voltage experiment, this method checks
      that the calculated voltage was within tolerance of the boundary
      condition. In addition, this method checks the anodic and cathodic
      Faradaic currents against the current at each time step and the
      solid-phase lithium conservation.

      If the verification returns ``False``, you can see which checks failed
      using ``plotflag``. Any subplots shaded grey failed their test.
      Failures generally suggest that the solver's relative and/or absolute
      tolerance should be adjusted, and/or that the discretization should be
      modified to adjust the mesh.

      :param plotflag: A flag to see plots showing the verification calculations. The
                       default is ``False``.
      :type plotflag: bool, optional
      :param rtol: Relative tolerance for comparisons. The default is 5e-3.
      :type rtol: float, optional
      :param atol: Absolute tolerance for comparisons. The default is 1e-3.
      :type atol: float, optional

      :returns: **checks** (*bool*) -- ``True`` is all checks are satisfied, ``False`` otherwise.



