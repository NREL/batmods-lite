:py:mod:`bmlite.P2D`
====================

.. py:module:: bmlite.P2D

.. autoapi-nested-parse::

   .. rubric:: Pseudo-2D Model Package

   A packaged pseudo-2D (P2D) model. Build a model using the ``Simulation`` class,
   and run any available experiments using its "run" methods, e.g. ``run_CC()``.
   The experiments return ``Solution`` class instances with post processing,
   plotting, and saving methods.



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

   dae/index.rst
   domains/index.rst
   postutils/index.rst
   roots/index.rst


Package Contents
----------------

Classes
~~~~~~~

.. autoapisummary::

   bmlite.P2D.Simulation



Functions
~~~~~~~~~

.. autoapisummary::

   bmlite.P2D.templates



.. py:class:: Simulation(yamlfile = 'graphite_nmc532')




   
   Make a P2D simulation capable of running various experiments.

   The initialization will add all of the battery attributes from the
   ``.yaml`` file under its ``bat``, ``el``, ``an``, ``sep``, and ``ca``
   attributes. The ``pre()`` method runs at the end of the initialization
   to add dependent parameters, including the mesh, algebraic indices,
   etc. to the simulation instance. This only happens in ``__init__``,
   which has some implications if the user modifies parameters after
   initialization (see the warning below).

   :param yamlfile: An absolute or relative path to the ``.yaml`` file that defines the
                    battery properties. The ``.yaml`` extension will be added to the
                    end of the string if it is not already there. The default is
                    ``'default_P2D'``, which loads an internal file from the ``bmlite``
                    library.
   :type yamlfile: str, optional

   .. warning::

      The user may choose to modify parameters after loading in a ``.yaml``
      file, however, they will need to manually re-run the ``pre()`` method
      if they do so. Otherwise, the dependent parameters may not be
      consistent with the user-defined inputs.

   .. seealso::

      :obj:`bmlite.P2D.templates`
          Get help making your own ``.yaml`` file by starting with the default template.

   .. py:method:: copy()

      Create a copy of the Simulation instance.

      :returns: **sim** (*P2D Simulation object*) -- A unique copy (stored separately in memory) of the Simulation
                instance.


   .. py:method:: j_pattern()

      Plot the Jacobian pattern.

      Runs the ``bmlite.P2D.dae.bandwidth`` function to determine and plot
      the Jacobian pattern.

      :returns: * **lband** (*int*) -- Lower bandwidth from the residual function's Jacobian pattern.
                * **uband** (*int*) -- Upper bandwidth from the residual function's Jacobian pattern.

      .. seealso:: :obj:`bmlite.P2D.dae.bandwidth`


   .. py:method:: pre()

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


   .. py:method:: run_CC(exp, **kwargs)

      Runs a constant current experiment specified by the details given in
      the experiment dictionary ``exp``.

      :param exp: The constant current experimental details. Required keys and
                  descriptions are listed below:

                  =========== ==============================================
                  Key         Value [units] (*type*)
                  =========== ==============================================
                  C_rate      C-rate (+ charge, - discharge) [1/h] (*float*)
                  t_min       minimum time [s] (*float*)
                  t_max       maximum time [s] (*float*)
                  Nt          number of time discretizations [-] (*int*)
                  =========== ==============================================
      :type exp: dict
      :param \*\*kwargs: The keyword arguments specify the Sundials IDA solver options. A
                         partial list of options/defaults is given below:

                         ============= =================================================
                         Key           Description (*type* or {options}, default)
                         ============= =================================================
                         rtol          relative tolerance (*float*, 1e-6)
                         atol          absolute tolerance (*float*, 1e-9)
                         linsolver     linear solver (``{'dense', 'band'}``, ``'band'``)
                         lband         width of the lower band (*int*, ``self.lband``)
                         uband         width of the upper band (*int*, ``self.uband``)
                         max_t_step    maximum time step (*float*, 0. -> unrestricted)
                         rootfn        root/event function (*Callable*, ``None``)
                         nr_rootfns    number of events in ``'rootfn'`` (*int*, 0)
                         ============= =================================================
      :type \*\*kwargs: dict, optional

      :returns: **sol** (*CCSolution object*) -- Solution class with the returned variable values, messages, exit
                flags, etc. from the IDA solver. The returned ``CCSolution``
                instance includes post processing, plotting, and saving methods.

      .. seealso:: :obj:`bmlite.IDASolver`, :obj:`bmlite.P2D.solutions.CCSolution`


   .. py:method:: run_CP(exp, **kwargs)

      Runs a constant power experiment specified by the details given in
      the experiment dictionary ``exp``.

      :param exp: The constant power experimental details. Required keys and
                  descriptions are listed below:

                  ======== ========================================================
                  Key      Value [units] (*type*)
                  ======== ========================================================
                  P_ext    external power (+ charge, - discharge) [W/m^2] (*float*)
                  t_min    minimum time [s] (*float*)
                  t_max    maximum time [s] (*float*)
                  Nt       number of time discretizations [-] (*int*)
                  ======== ========================================================
      :type exp: dict
      :param \*\*kwargs: The keyword arguments specify the Sundials IDA solver options. A
                         partial list of options/defaults is given below:

                         ============= =================================================
                         Key           Description (*type* or {options}, default)
                         ============= =================================================
                         rtol          relative tolerance (*float*, 1e-6)
                         atol          absolute tolerance (*float*, 1e-9)
                         linsolver     linear solver (``{'dense', 'band'}``, ``'band'``)
                         lband         width of the lower band (*int*, ``self.lband``)
                         uband         width of the upper band (*int*, ``self.uband``)
                         max_t_step    maximum time step (*float*, 0. -> unrestricted)
                         rootfn        root/event function (*Callable*, ``None``)
                         nr_rootfns    number of events in ``'rootfn'`` (*int*, 0)
                         ============= =================================================
      :type \*\*kwargs: dict, optional

      :returns: **sol** (*CPSolution object*) -- Solution class with the returned variable values, messages, exit
                flags, etc. from the IDA solver. The returned ``CPSolution``
                instance includes post processing, plotting, and saving methods.

      .. seealso:: :obj:`bmlite.IDASolver`, :obj:`bmlite.P2D.solutions.CPSolution`


   .. py:method:: run_CV(exp, **kwargs)

      Runs a constant voltage experiment specified by the details given in
      the experiment dictionary ``exp``.

      :param exp: The constant voltage experimental details. Required keys and
                  descriptions are listed below:

                  =========== ==========================================
                  Key         Value [units] (*type*)
                  =========== ==========================================
                  V_ext       externally applied voltage [V] (*float*)
                  t_min       minimum time [s] (*float*)
                  t_max       maximum time [s] (*float*)
                  Nt          number of time discretizations [-] (*int*)
                  =========== ==========================================
      :type exp: dict
      :param \*\*kwargs: The keyword arguments specify the Sundials IDA solver options. A
                         partial list of options/defaults is given below:

                         ============= =================================================
                         Key           Description (*type* or {options}, default)
                         ============= =================================================
                         rtol          relative tolerance (*float*, 1e-6)
                         atol          absolute tolerance (*float*, 1e-9)
                         linsolver     linear solver (``{'dense', 'band'}``, ``'band'``)
                         lband         width of the lower band (*int*, ``self.lband``)
                         uband         width of the upper band (*int*, ``self.uband``)
                         max_t_step    maximum time step (*float*, 0. -> unrestricted)
                         rootfn        root/event function (*Callable*, ``None``)
                         nr_rootfns    number of events in ``'rootfn'`` (*int*, 0)
                         ============= =================================================
      :type \*\*kwargs: dict, optional

      :returns: **sol** (*CVSolution object*) -- Solution class with the returned variable values, messages, exit
                flags, etc. from the IDA solver. The returned ``CVSolution``
                instance includes post processing, plotting, and saving methods.

      .. seealso:: :obj:`bmlite.IDASolver`, :obj:`bmlite.P2D.solutions.CVSolution`



.. py:function:: templates(sim = None, exp = None)

   Print simulation and/or experiment templates. If both ``sim`` and ``exp``
   are ``None``, a list of available templates will be printed. Otherwise, if
   a name or index is given, that template will print to the console.

   :param sim: Simulation template file name or index. The default is ``None``.
   :type sim: str | int, optional
   :param exp: Experiment template file name or index. The default is ``None``.
   :type exp: str | int, optional

   :returns: *None.*


