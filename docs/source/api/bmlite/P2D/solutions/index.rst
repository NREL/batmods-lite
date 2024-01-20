:py:mod:`bmlite.P2D.solutions`
==============================

.. py:module:: bmlite.P2D.solutions


Package Contents
----------------

Classes
~~~~~~~

.. autoapisummary::

   bmlite.P2D.solutions.BaseSolution
   bmlite.P2D.solutions.CCSolution




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

      :returns: **onroot** (*bool*) -- ``True`` if exit on root(event), ``False`` otherwise.

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

                  =========== =============================================
                  Key         Value (type)
                  =========== =============================================
                  t           solution times (1D array)
                  y           solution variables (2D array)
                  ydot        solution variable time derivatives (2D array)
                  success     overall solver exit status (bool)
                  onroot      onroot solver exit status (bool)
                  message     solver exit message (str)
                  solvetime   solver integration time, in seconds (float)
                  =========== =============================================
      :type sol: dict

      :returns: *None.*


   .. py:method:: ida_fill(sol: object, solvetime: float) -> None

      Fill the instance attributes using the SolverReturn object from the
      Sundials IDA solver.

      :param sol: Sundials IDA SolverReturn object with the following attributes:

                  =========== =============================================
                  Attribute   Value (type)
                  =========== =============================================
                  t           solution times (1D array)
                  y           solution variables (2D array)
                  ydot        solution variable time derivatives (2D array)
                  success     overall solver exit status (bool)
                  onroot      onroot solver exit status (bool)
                  message     solver exit message (str)
                  =========== =============================================
      :type sol: object
      :param solvetime: Solver integration time, in seconds.
      :type solvetime: float

      :returns: *None.*


   .. py:method:: report() -> None

      Prints the experiment detials and solution success report.

      :returns: *None.*


   .. py:method:: solvetime(units: str = 's') -> str

      Print solve time (not including pre/post processing).

      :param units: Units for printed time (``'s'``, ``'min'``, or ``'h'``). The
                    default is ``'s'``.
      :type units: str, optional

      :returns: **solvetime** (*str*) -- Time for Sundials IDA to solve problem in [units].


   .. py:method:: to_dict() -> dict

      Output a dictionary with key/value pairs corresponding to the instance
      attributes and values listed below.

      :returns: **sol** (*dict*) -- Solution dictionary with the following key/value pairs:

                =========== =============================================
                Key         Value (type)
                =========== =============================================
                t           solution times (1D array)
                y           solution variables (2D array)
                ydot        solution variable time derivatives (2D array)
                success     overall solver exit status (bool)
                onroot      onroot solver exit status (bool)
                message     solver exit message (str)
                solvetime   solver integration time, in seconds (float)
                =========== =============================================



.. py:class:: CCSolution(sim: object, exp: dict)




   Constant current solution for P2D simuations.

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

      Class name. Overwrites ``classname()`` from ``BaseSolution``.

      :returns: **classname** (*str*) -- Name of current class.

   .. py:method:: plot(*args) -> None


   .. py:method:: post() -> None


   .. py:method:: slice_and_save(savename: str, overwrite: bool = False) -> None

      Save a ``.npz`` file with all spatial, time, and state variables
      separated into 1D, 2D, and 3D arrays. The keys are given below.
      The index order of the 2D and 3D arrays is given with the value
      descriptions.

      ========= ==========================================================
      Key       Value [units] (type)
      ========= ==========================================================
      x_a       x mesh in anode [m] (1D array)
      x_s       x mesh in separator [m] (1D array)
      x_c       x mesh in cathode [m] (1D array)
      x         stacked x mesh for an, sep, and ca [m] (1D array)
      r_a       r mesh for anode particles [m] (1D array)
      r_c       r mesh for cathode particles [m] (1D array)
      t         saved solution times [s] (1D array)
      phie_a    electrolyte potentials at t, x_a [V] (2D array)
      phis_a    electrode potentials at t, x_a [V] (2D array)
      ce_a      electrolyte Li+ at t, x_a [kmol/m^3] (2D array)
      cs_a      electrode Li at t, x_a, r_a [kmol/m^3] (3D array)
      phie_s    electrolyte potentials at t, x_s [V] (2D array)
      ce_s      electrolyte Li+ at t, x_s [kmol/m^3] (2D array)
      phie_c    electrolyte potentials at t, x_c [V] (2D array)
      phis_c    electrode potentials at t, x_c [V] (2D array)
      ce_c      electrolyte Li+ at t, x_c [kmol/m^3] (2D array)
      cs_c      electrode Li at t, x_c, r_c [kmol/m^3] (3D array)
      phie      electrolyte potentials at t, x [V] (2D array)
      ce        electrolyte Li+ at t, x [kmol/m^3] (2D array)
      ie        electrolyte current at t, x boundarys [A/m^2] (2D array)
      j_a       anode Faradaic current at t, x_a [kmol/m^2/s] (2D array)
      j_c       cathode Faradaic current at t, x_c [kmol/m^2/s] (2D array)
      ========= ==========================================================

      :param savename: Either a file name or the absolute/relative file path. The ``.npz``
                       extension will be added to the end of the string if it is not
                       already there. If only the file name is given, the file will be
                       saved in the user's current working directory.
      :type savename: str
      :param overwrite: A flag to overwrite and existing ``.npz`` file with the same name
                        if one exists. The default is ``False``.
      :type overwrite: bool, optional

      :returns: *None.*



