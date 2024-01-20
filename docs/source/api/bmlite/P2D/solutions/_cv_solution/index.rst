:orphan:

:py:mod:`bmlite.P2D.solutions._cv_solution`
===========================================

.. py:module:: bmlite.P2D.solutions._cv_solution


Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   bmlite.P2D.solutions._cv_solution.CVSolution




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



