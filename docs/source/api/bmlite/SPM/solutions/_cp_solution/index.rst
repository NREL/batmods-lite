:orphan:

:py:mod:`bmlite.SPM.solutions._cp_solution`
===========================================

.. py:module:: bmlite.SPM.solutions._cp_solution


Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   bmlite.SPM.solutions._cp_solution.CPSolution




.. py:class:: CPSolution(sim: object, exp: dict)




   Constant power solution for SPM simuations.

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


   .. py:method:: slice_and_save(savename: str, overwrite: bool = False) -> None

      Save a ``.npz`` file with all spatial, time, and state variables
      separated into 1D and 2D arrays. The keys are given below. The index
      order of the 2D arrays is given with the value descriptions.

      ========= =====================================================
      Key       Value [units] (type)
      ========= =====================================================
      r_a       r mesh for anode particles [m] (1D array)
      r_c       r mesh for cathode particles [m] (1D array)
      t         saved solution times [s] (1D array)
      phis_a    anode electrode potentials at t [V] (1D array)
      cs_a      electrode Li at t, r_a [kmol/m^3] (2D array)
      phis_c    cathode electrode potentials at t [V] (1D array)
      cs_c      electrode Li at t, r_c [kmol/m^3] (2D array)
      phie      electrolyte potentials at t [V] (1D array)
      j_a       anode Faradaic current at t [kmol/m^2/s] (1D array)
      j_c       cathode Faradaic current at t [kmol/m^2/s] (1D array)
      ========= =====================================================

      :param savename: Either a file name or the absolute/relative file path. The ``.npz``
                       extension will be added to the end of the string if it is not
                       already there. If only the file name is given, the file will be
                       saved in the user's current working directory.
      :type savename: str
      :param overwrite: A flag to overwrite and existing ``.npz`` file with the same name
                        if one exists. The default is ``False``.
      :type overwrite: bool, optional

      :returns: *None.*



