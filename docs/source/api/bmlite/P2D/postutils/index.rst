bmlite.P2D.postutils
====================

.. py:module:: bmlite.P2D.postutils

.. autoapi-nested-parse::

   .. rubric:: Post-processing Utilities

   This module contains all post-processing functions for the P2D package. The
   available post-processing options for a given experiment are specific to that
   experiment. Therefore, not all ``Solution`` classes may have access to all of
   the following functions.



Functions
---------

.. autoapisummary::

   bmlite.P2D.postutils.electrolyte
   bmlite.P2D.postutils.intercalation
   bmlite.P2D.postutils.pixels
   bmlite.P2D.postutils.post
   bmlite.P2D.postutils.potentials


Module Contents
---------------

.. py:function:: electrolyte(soln)

   Plots electrolyte Li-ion concentration profiles vs. time.

   :param soln: A pseudo-2D model solution object.
   :type soln: P2D Solution object

   :returns: *None.*


.. py:function:: intercalation(soln)

   Plots anode and cathode particle intercalation profiles vs. time.

   :param soln: A pseudo-2D model solution object.
   :type soln: P2D Solution object

   :returns: *None.*


.. py:function:: pixels(soln)

   Makes pixel plots for most 2D (space/time) variables.

   :param soln: A pseudo-2D model solution object.
   :type soln: P2D Solution object

   :returns: *None.*


.. py:function:: post(soln)

   Run post processing to determine secondary variables.

   :param soln: A pseudo-2D model solution object.
   :type soln: Solution

   :returns: **postvars** (*dict*) -- Post processed variables, as described below.

             =========== ========================================================
             Key         Value [units] (*type*)
             =========== ========================================================
             div_i_an    divergence of current at t, x_a [A/m^3] (*2D array*)
             div_i_sep   divergence of current at t, x_s [A/m^3] (*2D array*)
             div_i_ca    divergence of current at t, x_c [A/m^3] (*2D array*)
             sdot_an     Faradaic current at t, x_a [kmol/m^2/s] (*1D array*)
             sdot_ca     Faradaic current at t, x_c [kmol/m^2/s] (*1D array*)
             sum_ip      ``i_ed + i_el`` at t, xp interfaces [A/m^2] (*2D array*)
             i_el_x      ``i_el`` at t, x interfaces [A/m^2] (*2D array*)
             =========== ========================================================

   .. seealso:: :obj:`~bmlite.P2D.solutions.StepSolution`, :obj:`~bmlite.P2D.solutions.CycleSolution`


.. py:function:: potentials(soln)

   Plots anode, electrolyte, and cathode potentials vs. time and space.

   :param soln: A pseudo-2D model solution object.
   :type soln: P2D Solution object

   :returns: *None.*


