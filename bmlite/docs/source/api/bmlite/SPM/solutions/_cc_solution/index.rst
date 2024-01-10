:orphan:

:py:mod:`bmlite.SPM.solutions._cc_solution`
===========================================

.. py:module:: bmlite.SPM.solutions._cc_solution


Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   bmlite.SPM.solutions._cc_solution.CCSolution




.. py:class:: CCSolution(sim: object, exp: dict)




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

      Class name. Overwrites ``classname()`` from ``BaseSolution``.

      :returns: **classname** (*str*) -- Name of current class.

   .. py:method:: plot(*args) -> None


   .. py:method:: post() -> None



