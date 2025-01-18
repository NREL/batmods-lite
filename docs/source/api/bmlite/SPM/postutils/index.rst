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
   bmlite.SPM.postutils.potentials


Module Contents
---------------

.. py:function:: intercalation(sol)

   Plots anode and cathode particle intercalation profiles vs. time.

   :param sol: A single particle model solution object.
   :type sol: SPM Solution object

   :returns: *None.*


.. py:function:: pixels(sol)

   Makes pixel plots for most 2D (space/time) variables.

   :param sol: A single particle model solution object.
   :type sol: SPM Solution object

   :returns: *None.*


.. py:function:: potentials(sol)

   Plots anode, electrolyte, and cathode potentials vs. time.

   :param sol: A single particle model solution object.
   :type sol: SPM Solution object

   :returns: *None.*


