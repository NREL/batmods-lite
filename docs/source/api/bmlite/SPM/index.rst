:py:mod:`bmlite.SPM`
====================

.. py:module:: bmlite.SPM

.. autoapi-nested-parse::

   .. rubric:: Single Particle Model Package

   A packaged single particle model (SPM). Build a model using the ``Simulation``
   class, and run any available experiments using its "run" methods, e.g.
   ``run_CC()``. The experiments return ``Solution`` class instances with post
   processing, plotting, and saving methods.



Subpackages
-----------
.. toctree::
   :titlesonly:
   :maxdepth: 3

   solutions/index.rst


Submodules
----------
.. toctree::
   :titlesonly:
   :maxdepth: 1

   builder/index.rst
   dae/index.rst
   postutils/index.rst


Package Contents
----------------

Classes
~~~~~~~

.. autoapisummary::

   bmlite.SPM.Simulation



Functions
~~~~~~~~~

.. autoapisummary::

   bmlite.SPM.load
   bmlite.SPM.templates



.. py:class:: Simulation(yamlfile: str = 'default_SPM')




   
   Make a SPM simulation capable of running various experiments.

   The initialization will add all of the battery attributes from the
   ``.yaml`` file under its ``bat``, ``el``, ``an``, and ``ca``
   attributes. The ``pre()`` method runs at the end of the initialization
   to add dependent parameters, including the mesh, algebraic indices,
   etc. to the simulation instance. This only happens in ``__init__``,
   which has some implications if the user modifies parameters after
   initialization (see the warning below).

   :param yamlfile: An absolute or relative path to the ``.yaml`` file that defines the
                    battery properties. The ``.yaml`` extension will be added to the
                    end of the string if it is not already there. The default is
                    ``'default_SPM'``, which loads an internal file from the ``bmlite``
                    library.
   :type yamlfile: str, optional

   .. warning::

      The user may choose to modify parameters after loading in a ``.yaml``
      file, however, they will need to manually re-run the ``pre()`` method
      if they do so. Otherwise, the dependent parameters may not be
      consistent with the user-defined inputs.

   .. rubric:: Notes

   * The ``'default_SPM.yaml'`` file can be accessed via the
     ``templates()`` method. Please reference this for help building your
     own ``.yaml`` file.

   .. py:method:: j_pattern() -> None

      Plot the Jacobian pattern.

      Runs the ``bmlite.SPM.dae.bandwidth`` function to determine and plot
      the Jacobian pattern.

      :returns: *None.*

      .. seealso:: :obj:`bmlite.SPM.dae.bandwidth`


   .. py:method:: pre() -> None

      Pre-process the dependent parameters.

      The dependent parameters include ``A_s``, ``eps_s``, ``eps_AM``,
      ``sigma_s``, and setting the ``material`` classes for each domain. In
      addition, this method determines the mesh, pointers, algebraic indices,
      bandwidth, and initial solution. ``pre()`` is automatically executed
      in the ``__init__()`` method which has some implications if the user
      modifies parameters after initialization (see the warning below).

      :returns: *None.*

      .. warning::

         The user may choose to modify parameters after loading in a ``.yaml``
         file, however, they will need to manually re-run the ``pre()`` method
         if they do so. Otherwise, the dependent parameters may not be
         consistent with the user-defined inputs.


   .. py:method:: run_CC(exp: dict, **kwargs) -> object

      Runs a constant current experiment specified by the details given in
      the experiment dictionary ``exp``.

      :param exp: The constant current experimental details. Required keys and
                  descriptions are listed below:

                  =========== ==============================================
                  Key         Value [units] (type)
                  =========== ==============================================
                  C_rate      C-rate (+ charge, - discharge) [1/h] (*float*)
                  t_min       minimum time [s] (*float*)
                  t_max       maximum time [s] (*float*)
                  Nt          number of time discretizations [-] (*int*)
                  =========== ==============================================
      :type exp: dict
      :param \*\*kwargs: The keyword arguments specify the Sundials IDA solver options. A
                         partial list of options/defaults is given below:

                         =============== =================================================
                         Key             Description (type or options, default)
                         =============== =================================================
                         rtol            relative tolerance (*float*, 1e-6)
                         atol            absolute tolerance (*float*, 1e-9)
                         linsolver       linear solver (``{'dense', 'band'}``, ``'band'``)
                         lband           width of the lower band (*int*, ``self.lband``)
                         uband           width of the upper band (*int*, ``self.uband``)
                         max_step_size   maximum time step (*float*, 0. -> unrestricted)
                         rootfn          root/event function (*Callable*, ``None``)
                         nr_rootfns      number of events in ``'rootfn'`` (*int*, 0)
                         =============== =================================================
      :type \*\*kwargs: dict, optional

      :returns: **sol** (*CCSolution object*) -- Solution class with the returned variable values, messages, exit
                flags, etc. from the IDA solver. The returned ``CCSolution``
                instance includes post processing, plotting, and saving methods.

      .. seealso:: :obj:`bmlite.IDASolver`, :obj:`bmlite.SPM.solutions.CCSolution`


   .. py:method:: run_CP(exp: dict, **kwargs) -> object

      Runs a constant power experiment specified by the details given in
      the experiment dictionary ``exp``.

      :param exp: The constant power experimental details. Required keys and
                  descriptions are listed below:

                  ======== ========================================================
                  Key      Value [units] (type)
                  ======== ========================================================
                  P_ext    external power (+ charge, - discharge) [W/m^2] (*float*)
                  t_min    minimum time [s] (*float*)
                  t_max    maximum time [s] (*float*)
                  Nt       number of time discretizations [-] (*int*)
                  ======== ========================================================
      :type exp: dict
      :param \*\*kwargs: The keyword arguments specify the Sundials IDA solver options. A
                         partial list of options/defaults is given below:

                         =============== =================================================
                         Key             Description (type or options, default)
                         =============== =================================================
                         rtol            relative tolerance (*float*, 1e-6)
                         atol            absolute tolerance (*float*, 1e-9)
                         linsolver       linear solver (``{'dense', 'band'}``, ``'band'``)
                         lband           width of the lower band (*int*, ``self.lband``)
                         uband           width of the upper band (*int*, ``self.uband``)
                         max_step_size   maximum time step (*float*, 0. -> unrestricted)
                         rootfn          root/event function (*Callable*, ``None``)
                         nr_rootfns      number of events in ``'rootfn'`` (*int*, 0)
                         =============== =================================================
      :type \*\*kwargs: dict, optional

      :returns: **sol** (*CPSolution object*) -- Solution class with the returned variable values, messages, exit
                flags, etc. from the IDA solver. The returned ``CPSolution``
                instance includes post processing, plotting, and saving methods.

      .. seealso:: :obj:`bmlite.IDASolver`, :obj:`bmlite.SPM.solutions.CPSolution`


   .. py:method:: run_CV(exp: dict, **kwargs) -> object

      Runs a constant voltage experiment specified by the details given in
      the experiment dictionary ``exp``.

      :param exp: The constant voltage experimental details. Required keys and
                  descriptions are listed below:

                  =========== ==========================================
                  Key         Value [units] (type)
                  =========== ==========================================
                  V_ext       externally applied voltage [V] (*float*)
                  t_min       minimum time [s] (*float*)
                  t_max       maximum time [s] (*float*)
                  Nt          number of time discretizations [-] (*int*)
                  =========== ==========================================
      :type exp: dict
      :param \*\*kwargs: The keyword arguments specify the Sundials IDA solver options. A
                         partial list of options/defaults is given below:

                         =============== =================================================
                         Key             Description (type or options, default)
                         =============== =================================================
                         rtol            relative tolerance (*float*, 1e-6)
                         atol            absolute tolerance (*float*, 1e-9)
                         linsolver       linear solver (``{'dense', 'band'}``, ``'band'``)
                         lband           width of the lower band (*int*, ``self.lband``)
                         uband           width of the upper band (*int*, ``self.uband``)
                         max_step_size   maximum time step (*float*, 0. -> unrestricted)
                         rootfn          root/event function (*Callable*, ``None``)
                         nr_rootfns      number of events in ``'rootfn'`` (*int*, 0)
                         =============== =================================================
      :type \*\*kwargs: dict, optional

      :returns: **sol** (*CVSolution object*) -- Solution class with the returned variable values, messages, exit
                flags, etc. from the IDA solver. The returned ``CVSolution``
                instance includes post processing, plotting, and saving methods.

      .. seealso:: :obj:`bmlite.IDASolver`, :obj:`bmlite.SPM.solutions.CVSolution`



.. py:function:: load(loadname: str) -> tuple[object]

   Load a previous SPM simulation and solution pair.

   The ``Solution`` classes have a save method that create a directory
   with saved ``.yaml`` and ``.npz`` files which can be used to reconstruct
   the ``sim`` and ``sol`` objects from a previous experiment.

   :param loadname: The absolute or relative path to a directory with previously saved
                    solution files.
   :type loadname: str

   :returns: * **sim** (*SPM Simulation object*) -- An initialized SPM simulation instance.
             * **sol** (*SPM Solution object*) -- An initialized and filled SPM solution instance. The solution class is
               determined from information in the previously saved ``.yaml`` files.

   .. seealso:: :obj:`bmlite.SPM.Simulation`, :obj:`bmlite.SPM.solutions.BaseSolution.save`


.. py:function:: templates(sim: str | int = None, exp: str | int = None) -> None

   Print simulation and/or experiment templates. If both ``sim`` and ``exp``
   are ``None``, a list of available templates will be printed. Otherwise, if
   a name or index is given, that template will print to the console.

   :param sim: Simulation template file name or index. The default is ``None``.
   :type sim: str | int, optional
   :param exp: Experiment template file name or index. The default is ``None``.
   :type exp: str | int, optional

   :returns: *None.*


