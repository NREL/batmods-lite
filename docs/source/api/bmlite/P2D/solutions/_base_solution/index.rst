:orphan:

:py:mod:`bmlite.P2D.solutions._base_solution`
=============================================

.. py:module:: bmlite.P2D.solutions._base_solution


Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   bmlite.P2D.solutions._base_solution.BaseSolution




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



