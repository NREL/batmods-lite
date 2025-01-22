bmlite.SPM.postutils
====================

.. py:module:: bmlite.SPM.postutils

.. autoapi-nested-parse::

   .. rubric:: Post-processing Utilities

   This module contains all post-processing functions for the SPM package. The
   available post-processing options for a given experiment are specific to that
   experiment. Therefore, not all ``Solution`` classes may have access to all of
   the following functions.



Functions
---------

.. autoapisummary::

   bmlite.SPM.postutils.intercalation
   bmlite.SPM.postutils.pixels
   bmlite.SPM.postutils.post
   bmlite.SPM.postutils.potentials


Module Contents
---------------

.. py:function:: intercalation(soln)

   Plots anode and cathode particle intercalation profiles vs. time.

   :param soln: A single particle model solution object.
   :type soln: Solution

   :returns: *None.*

   .. seealso:: :obj:`~bmlite.SPM.solutions.StepSolution`, :obj:`~bmlite.SPM.solutions.CycleSolution`


.. py:function:: pixels(soln)

   Makes pixel plots for most 2D (space/time) variables.

   :param soln: A single particle model solution object.
   :type soln: SPM Solution object

   :returns: *None.*

   .. seealso:: :obj:`~bmlite.SPM.solutions.StepSolution`, :obj:`~bmlite.SPM.solutions.CycleSolution`


.. py:function:: post(soln)

   Run post processing to determine secondary variables.

   :param soln: A single particle model solution object.
   :type soln: Solution

   :returns: **postvars** (*dict*) -- Post processed variables, as described below.

             ========= ========================================================
             Key       Value [units] (*type*)
             ========= ========================================================
             sdot_an   anode Faradaic current at t [kmol/m^2/s] (*1D array*)
             sdot_ca   cathode Faradaic current at t [kmol/m^2/s] (*1D array*)
             ========= ========================================================

   .. seealso:: :obj:`~bmlite.SPM.solutions.StepSolution`, :obj:`~bmlite.SPM.solutions.CycleSolution`


.. py:function:: potentials(soln)

   Plots anode, electrolyte, and cathode potentials vs. time.

   :param soln: A single particle model solution object.
   :type soln: Solution

   :returns: *None.*

   .. seealso:: :obj:`~bmlite.SPM.solutions.StepSolution`, :obj:`~bmlite.SPM.solutions.CycleSolution`


