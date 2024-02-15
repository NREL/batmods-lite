:orphan:

:py:mod:`bmlite.P2D.solutions._cp_solution`
===========================================

.. py:module:: bmlite.P2D.solutions._cp_solution


Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   bmlite.P2D.solutions._cp_solution.CPSolution




.. py:class:: CPSolution(sim, exp)




   Constant power solution for P2D simulations.

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

   .. py:method:: verify(plotflag = False, rtol = 0.005, atol = 0.001)

      Verifies the solution is mathematically consistent.

      Specifically, for a constant power experiment, this method checks
      that the calculated power was within tolerance of the boundary
      condition. In addition, this method checks the anodic and cathodic
      Faradaic currents against the current at each time step, the
      through-current in each control volume, and the Li-ion and solid-phase
      lithium conservation.

      If the verification returns ``False``, you can see which checks failed
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



