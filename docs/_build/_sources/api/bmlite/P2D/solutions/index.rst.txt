:py:mod:`bmlite.P2D.solutions`
==============================

.. py:module:: bmlite.P2D.solutions

.. autoapi-nested-parse::

   .. rubric:: Solution Classes

   This package contains class definitions for each experimental solution. Most
   solutions inherit from :class:`~bmlite.P2D.solutions.BaseSolution`, however,
   their verification processes are unique.



Package Contents
----------------

Classes
~~~~~~~

.. autoapisummary::

   bmlite.P2D.solutions.BaseSolution
   bmlite.P2D.solutions.CCSolution
   bmlite.P2D.solutions.CPSolution
   bmlite.P2D.solutions.CVSolution




.. py:class:: BaseSolution(sim: object, exp: dict)




   Base methods for all P2D Solution classes.

   When initialized, a copy of the P2D simulation object and experimental
   details for the solution instance are stored, all other instance
   attributes are set to ``None``. These are filled later with a "fill"
   method (e.g., ``ida_fill()`` or ``dict_fill()``).

   :param sim: The P2D Simulation instance used to produce the solution.
   :type sim: P2D Simulation object
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

   .. py:method:: dict_fill(sol: dict) -> None

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


   .. py:method:: ida_fill(sol: object, solvetime: float) -> None

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


   .. py:method:: plot(*args: str) -> None

      Generates requested plots based on ``*args``.

      :param \*args: Use any number of the following arguments to see the described
                     plots:

                     ================= ===============================================
                     arg               Plot description
                     ================= ===============================================
                     'current'         current density [A/m^2] vs. time [s]
                     'voltage'         cell voltage [V] vs. time [s]
                     'power'           power density [W/m^2] vs. time [s]
                     'ivp'             combined current, voltage, and power plot
                     'potentials'      anode, cathode, and electrolyte potentials [V]
                     'electrolyte'     Li-ion concentration [kmol/m^3] vs. x and t
                     'intercalation'   anode/cathode particle Li fractions vs. r and t
                     'pixels'          pixel plots for most 2D variables
                     ================= ===============================================
      :type \*args: str

      :returns: *None.*


   .. py:method:: post() -> None


   .. py:method:: report() -> None

      Prints the experiment details and solution success report.

      :returns: *None.*


   .. py:method:: slice_and_save(savename: str, overwrite: bool = False) -> None

      Save a ``.npz`` file with all spatial, time, and state variables
      separated into 1D, 2D, and 3D arrays. The keys are given below.
      The index order of the 2D and 3D arrays is given with the value
      descriptions.

      ========= ====================================================
      Key       Value [units] (*type*)
      ========= ====================================================
      x_a       x mesh in anode [m] (*1D array*)
      x_s       x mesh in separator [m] (*1D array*)
      x_c       x mesh in cathode [m] (*1D array*)
      x         stacked x mesh for an, sep, and ca [m] (*1D array*)
      r_a       r mesh for anode particles [m] (*1D array*)
      r_c       r mesh for cathode particles [m] (*1D array*)
      t         saved solution times [s] (*1D array*)
      phie_a    electrolyte potentials at t, x_a [V] (*2D array*)
      phis_a    electrode potentials at t, x_a [V] (*2D array*)
      ce_a      electrolyte Li+ at t, x_a [kmol/m^3] (*2D array*)
      cs_a      electrode Li at t, x_a, r_a [kmol/m^3] (*3D array*)
      phie_s    electrolyte potentials at t, x_s [V] (*2D array*)
      ce_s      electrolyte Li+ at t, x_s [kmol/m^3] (*2D array*)
      phie_c    electrolyte potentials at t, x_c [V] (*2D array*)
      phis_c    electrode potentials at t, x_c [V] (*2D array*)
      ce_c      electrolyte Li+ at t, x_c [kmol/m^3] (*2D array*)
      cs_c      electrode Li at t, x_c, r_c [kmol/m^3] (*3D array*)
      phie      electrolyte potentials at t, x [V] (*2D array*)
      ce        electrolyte Li+ at t, x [kmol/m^3] (*2D array*)
      ie        ``i_el`` at t, x boundarys [A/m^2] (*2D array*)
      j_a       Faradaic current at t, x_a [kmol/m^2/s] (*2D array*)
      j_c       Faradaic current at t, x_c [kmol/m^2/s] (*2D array*)
      ========= ====================================================

      :param savename: Either a file name or the absolute/relative file path. The ``.npz``
                       extension will be added to the end of the string if it is not
                       already there. If only the file name is given, the file will be
                       saved in the user's current working directory.
      :type savename: str
      :param overwrite: A flag to overwrite an existing ``.npz`` file with the same name
                        if one exists. The default is ``False``.
      :type overwrite: bool, optional

      :returns: *None.*


   .. py:method:: solvetime(units: str = 's') -> str

      Print solve time (not including pre/post processing).

      :param units: Units for printed time ``{'s', 'min', 'h'}``. The default is
                    ``'s'``.
      :type units: str, optional

      :returns: **solvetime** (*str*) -- Time for Sundials IDA to solve problem in [units].


   .. py:method:: to_dict() -> dict

      Output a dictionary with key/value pairs corresponding to the instance
      attributes and values listed below.

      :returns: **sol** (*dict*) -- Solution dictionary with the following key/value pairs:

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


   .. py:method:: to_sliced_dict() -> dict

      Creates a dict with all spatial, time, and state variables
      separated into 1D, 2D, and 3D arrays. The keys are given below.
      The index order of the 2D and 3D arrays is given with the value
      descriptions.

      ========= ====================================================
      Key       Value [units] (*type*)
      ========= ====================================================
      x_a       x mesh in anode [m] (*1D array*)
      x_s       x mesh in separator [m] (*1D array*)
      x_c       x mesh in cathode [m] (*1D array*)
      x         stacked x mesh for an, sep, and ca [m] (*1D array*)
      r_a       r mesh for anode particles [m] (*1D array*)
      r_c       r mesh for cathode particles [m] (*1D array*)
      t         saved solution times [s] (*1D array*)
      phie_a    electrolyte potentials at t, x_a [V] (*2D array*)
      phis_a    electrode potentials at t, x_a [V] (*2D array*)
      ce_a      electrolyte Li+ at t, x_a [kmol/m^3] (*2D array*)
      cs_a      electrode Li at t, x_a, r_a [kmol/m^3] (*3D array*)
      phie_s    electrolyte potentials at t, x_s [V] (*2D array*)
      ce_s      electrolyte Li+ at t, x_s [kmol/m^3] (*2D array*)
      phie_c    electrolyte potentials at t, x_c [V] (*2D array*)
      phis_c    electrode potentials at t, x_c [V] (*2D array*)
      ce_c      electrolyte Li+ at t, x_c [kmol/m^3] (*2D array*)
      cs_c      electrode Li at t, x_c, r_c [kmol/m^3] (*3D array*)
      phie      electrolyte potentials at t, x [V] (*2D array*)
      ce        electrolyte Li+ at t, x [kmol/m^3] (*2D array*)
      ie        ``i_el`` at t, x boundarys [A/m^2] (*2D array*)
      j_a       Faradaic current at t, x_a [kmol/m^2/s] (*2D array*)
      j_c       Faradaic current at t, x_c [kmol/m^2/s] (*2D array*)
      ========= ====================================================

      :param None.:

      :returns: **sol_dict** (*dict*) -- A dictionary containing the solution.



.. py:class:: CCSolution(sim: object, exp: dict)




   Constant current solution for P2D simulations.

   Base: :class:`~bmlite.P2D.solutions.BaseSolution`

   When initialized, a copy of the P2D simulation object and experimental
   details for the solution instance are stored, all other instance
   attributes are set to ``None``. These are filled later with a "fill"
   method (e.g., ``ida_fill()`` or ``dict_fill()``).

   :param sim: The P2D Simulation instance used to produce the solution.
   :type sim: P2D Simulation object
   :param exp: Experiment dictionary. Specific key/value pairs are dependent on
               the experiment that was run.
   :type exp: dict

   .. py:property:: classname
      :type: str

      Class name. Overwrites ``BaseSolution``.

      :returns: **classname** (*str*) -- Name of current class.

   .. py:method:: verify(plotflag: bool = False, rtol: float = 0.005, atol: float = 0.001) -> bool

      Verifies the solution is mathematically consistent.

      Specifically, for a constant current experiment, this method checks
      that the calculated current was within tolerance of the boundary
      condition. In addition, this method checks the anodic and cathodic
      Faradaic currents against the current at each time step, the
      through-current in each control volume, and the Li-ion and solid-phase
      lithium conservation.

      If the verificaion returns ``False``, you can see which checks failed
      using ``plotflag``. Any subplots shaded grey failed their test.
      Failures generally suggest that the solver's relative and/or absolute
      tolerance should be adjusted, and/or that the discretization should be
      modified to adjust the mesh. Note that the third row of the figure
      shows conservation of charge in each control volume in each domain.
      This row is not included in the checks and will never shade grey, but
      is useful in debugging.

      :param plotflag: A flag to see plots showing the verification calculations. The
                       default is ``False``.
      :type plotflag: bool, optional
      :param rtol: Relative tolerance for comparisons. The default is 5e-3.
      :type rtol: float, optional
      :param atol: Absolute tolerance for comparisons. The default is 1e-3.
      :type atol: float, optional

      :returns: **checks** (*bool*) -- ``True`` is all checks are satisfied, ``False`` otherwise.



.. py:class:: CPSolution(sim: object, exp: dict)




   Constant power solution for P2D simuations.

   Base: :class:`~bmlite.P2D.solutions.BaseSolution`

   When initialized, a copy of the P2D simulation object and experimental
   details for the solution instance are stored, all other instance
   attributes are set to ``None``. These are filled later with a "fill"
   method (e.g., ``ida_fill()`` or ``dict_fill()``).

   :param sim: The P2D Simulation instance used to produce the solution.
   :type sim: P2D Simulation object
   :param exp: Experiment dictionary. Specific key/value pairs are dependent on
               the experiment that was run.
   :type exp: dict

   .. py:property:: classname
      :type: str

      Class name. Overwrites ``BaseSolution``.

      :returns: **classname** (*str*) -- Name of current class.

   .. py:method:: verify(plotflag: bool = False, rtol: float = 0.005, atol: float = 0.001) -> bool

      Verifies the solution is mathematically consistent.

      Specifically, for a constant power experiment, this method checks
      that the calculated power was within tolerance of the boundary
      condition. In addition, this method checks the anodic and cathodic
      Faradaic currents against the current at each time step, the
      through-current in each control volume, and the Li-ion and solid-phase
      lithium conservation.

      If the verificaion returns ``False``, you can see which checks failed
      using ``plotflag``. Any subplots shaded grey failed their test.
      Failures generally suggest that the solver's relative and/or absolute
      tolerance should be adjusted, and/or that the discretization should be
      modified to adjust the mesh. Note that the third row of the figure
      shows conservation of charge in each control volume in each domain.
      This row is not included in the checks and will never shade grey, but
      is useful in debugging.

      :param plotflag: A flag to see plots showing the verification calculations. The
                       default is ``False``.
      :type plotflag: bool, optional
      :param rtol: Relative tolerance for comparisons. The default is 5e-3.
      :type rtol: float, optional
      :param atol: Absolute tolerance for comparisons. The default is 1e-3.
      :type atol: float, optional

      :returns: **checks** (*bool*) -- ``True`` is all checks are satisfied, ``False`` otherwise.



.. py:class:: CVSolution(sim: object, exp: dict)




   Constant voltage solution for P2D simuations.

   Base: :class:`~bmlite.P2D.solutions.BaseSolution`

   When initialized, a copy of the P2D simulation object and experimental
   details for the solution instance are stored, all other instance
   attributes are set to ``None``. These are filled later with a "fill"
   method (e.g., ``ida_fill()`` or ``dict_fill()``).

   :param sim: The P2D Simulation instance used to produce the solution.
   :type sim: P2D Simulation object
   :param exp: Experiment dictionary. Specific key/value pairs are dependent on
               the experiment that was run.
   :type exp: dict

   .. py:property:: classname
      :type: str

      Class name. Overwrites ``BaseSolution``.

      :returns: **classname** (*str*) -- Name of current class.

   .. py:method:: verify(plotflag: bool = False, rtol: float = 0.005, atol: float = 0.001) -> bool

      Verifies the solution is mathematically consistent.

      Specifically, for a constant voltage experiment, this method checks
      that the calculated voltage was within tolerance of the boundary
      condition. In addition, this method checks the anodic and cathodic
      Faradaic currents against the current at each time step, the
      through-current in each control volume, and the Li-ion and solid-phase
      lithium conservation.

      If the verificaion returns ``False``, you can see which checks failed
      using ``plotflag``. Any subplots shaded grey failed their test.
      Failures generally suggest that the solver's relative and/or absolute
      tolerance should be adjusted, and/or that the discretization should be
      modified to adjust the mesh. Note that the third row of the figure
      shows conservation of charge in each control volume in each domain.
      This row is not included in the checks and will never shade grey, but
      is useful in debugging.

      :param plotflag: A flag to see plots showing the verification calculations. The
                       default is ``False``.
      :type plotflag: bool, optional
      :param rtol: Relative tolerance for comparisons. The default is 5e-3.
      :type rtol: float, optional
      :param atol: Absolute tolerance for comparisons. The default is 1e-3.
      :type atol: float, optional

      :returns: **checks** (*bool*) -- ``True`` is all checks are satisfied, ``False`` otherwise.



