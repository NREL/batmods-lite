bmlite.SPM
==========

.. py:module:: bmlite.SPM

.. autoapi-nested-parse::

   .. rubric:: Single Particle Model Package

   A packaged single particle model (SPM). Build a model using the ``Simulation``
   class, and run any available experiments using its "run" methods, e.g.
   ``run_CC()``. The experiments return ``Solution`` class instances with post
   processing, plotting, and saving methods.



Submodules
----------

.. toctree::
   :maxdepth: 1

   /api/bmlite/SPM/dae/index
   /api/bmlite/SPM/domains/index
   /api/bmlite/SPM/postutils/index


Classes
-------

.. autoapisummary::

   bmlite.SPM.CycleSolution
   bmlite.SPM.Simulation
   bmlite.SPM.StepSolution


Functions
---------

.. autoapisummary::

   bmlite.SPM.templates


Package Contents
----------------

.. py:class:: CycleSolution(*soln, t_shift = 0.001)



   All-step solution.

   A solution instance with all experiment steps stitch together into
   a single cycle.

   :param \*soln: All unpacked StepSolution instances to stitch together. The given
                  steps should be given in the same sequential order that they were
                  run.
   :type \*soln: StepSolution
   :param t_shift: Time (in seconds) to shift step solutions by when stitching them
                   together. If zero the end time of each step overlaps the starting
                   time of its following step. The default is 1e-3.
   :type t_shift: float


   .. py:method:: get_steps(idx)

      Return a subset of the solution.

      :param idx: The step index (int) or first/last indices (tuple) to return.
      :type idx: int | tuple

      :returns: :class:`StepSolution` | :class:`CycleSolution` -- The returned solution subset. A StepSolution is returned if 'idx'
                is an int, and a CycleSolution will be returned for the range of
                requested steps when 'idx' is a tuple.



   .. py:property:: solvetime
      :type: str

      Print a statement specifying how long IDASolver spent integrating.

      :returns: **solvetime** (*str*) -- An f-string with the total solver integration time in seconds.


.. py:class:: Simulation(yamlfile = 'graphite_nmc532')



   
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

   .. seealso::

      :obj:`bmlite.SPM.templates`
          Get help making your own ``.yaml`` file by starting with the default template.


   .. py:method:: copy()

      Create a copy of the Simulation instance.

      :returns: **sim** (*SPM Simulation object*) -- A unique copy (stored separately in memory) of the Simulation
                instance.



   .. py:method:: j_pattern()

      Plot the Jacobian pattern.

      Runs the ``bmlite.SPM.dae.bandwidth`` function to determine and plot
      the Jacobian pattern.

      :returns: * **lband** (*int*) -- Lower bandwidth from the residual function's Jacobian pattern.
                * **uband** (*int*) -- Upper bandwidth from the residual function's Jacobian pattern.

      .. seealso:: :obj:`bmlite.SPM.dae.bandwidth`



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



   .. py:method:: run(expr, reset_state = True, t_shift = 0.001)

      Run a full experiment.

      :param expr: An experiment instance.
      :type expr: Experiment
      :param reset_state: If True (default), the internal state of the model will be reset
                          back to a rested condition at 'soc0' at the end of all steps. When
                          False, the state does not reset. Instead it will update to match
                          the final state of the last experimental step.
      :type reset_state: bool
      :param t_shift: Time (in seconds) to shift step solutions by when stitching them
                      together. If zero the end time of each step overlaps the starting
                      time of its following step. The default is 1e-3.
      :type t_shift: float

      :returns: :class:`~bmlite.SPM.CycleSolution` -- A stitched solution with all experimental steps.

      .. warning::

         The default behavior resets the model's internal state back to a rested
         condition at 'soc0' by calling the ``pre()`` method at the end of all
         steps. This means that if you run a second experiment afterward, it
         will not start where the previous one left off. Instead, it will start
         from the original rested condition that the model initialized with. You
         can bypass this by using ``reset_state=False``, which keeps the state
         at the end of the final experimental step.

      .. seealso::

         :obj:`Experiment`
             Build an experiment.

         :obj:`CycleSolution`
             Wrapper for an all-steps solution.



   .. py:method:: run_step(expr, stepidx)

      Run a single experimental step.

      :param expr: An experiment instance.
      :type expr: Experiment
      :param stepidx: Step index to run. The first step has index 0.
      :type stepidx: int

      :returns: :class:`~bmlite.SPM.StepSolution` -- Solution to the experimental step.

      .. warning::

         The model's internal state is changed at the end of each experimental
         step. Consequently, you should not run steps out of order. You should
         always start with ``stepidx = 0`` and then progress to the subsequent
         steps afterward. Run ``pre()`` after your last step to reset the state
         back to a rested condition at 'soc0', if needed. Alternatively, you
         can continue running experiments back-to-back without a pre-processing
         in between if you want the following experiment to pick up from the
         same state that the last experiment ended.

      .. seealso::

         :obj:`Experiment`
             Build an experiment.

         :obj:`StepSolution`
             Wrapper for a single-step solution.

      .. rubric:: Notes

      Using the ``run()`` loops through all steps in an experiment and then
      stitches their solutions together. Most of the time, this is more
      convenient. However, advantages for running step-by-step is that it
      makes it easier to fine tune solver options, and allows for analyses
      or control decisions in the middle of an experiment.



.. py:class:: StepSolution(sim, idasoln, timer)



   Single-step solution.

   A solution instance for a single experimental step.

   :param sim: The simulation instance that was run to produce the solution.
   :type sim: Simulation
   :param idasoln: The unformatted solution returned by IDASolver.
   :type idasoln: IDAResult
   :param timer: Amount of time it took for IDASolver to perform the integration.
   :type timer: float


   .. py:property:: solvetime
      :type: str

      Print a statement specifying how long IDASolver spent integrating.

      :returns: **solvetime** (*str*) -- An f-string with the solver integration time in seconds.


.. py:function:: templates(sim = None)

   Print simulation templates. If ``sim`` is ``None``, a list of available
   templates will be printed. Otherwise, if a name or index is given, that
   template will print to the console.

   :param sim: Simulation template file name or index. The default is ``None``.
   :type sim: str | int, optional

   :returns: *None.*


