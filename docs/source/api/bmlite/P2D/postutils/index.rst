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
   bmlite.P2D.postutils.potentials


Module Contents
---------------

.. py:function:: electrolyte(sol)

   Plots electrolyte Li-ion concentration profiles vs. time.

   :param sol: A pseudo-2D model solution object.
   :type sol: P2D Solution object

   :returns: *None.*


.. py:function:: intercalation(sol)

   Plots anode and cathode particle intercalation profiles vs. time.

   :param sol: A pseudo-2D model solution object.
   :type sol: P2D Solution object

   :returns: *None.*


.. py:function:: pixels(sol)

   Makes pixel plots for most 2D (space/time) variables.

   :param sol: A pseudo-2D model solution object.
   :type sol: P2D Solution object

   :returns: *None.*


.. py:function:: potentials(sol)

   Plots anode, electrolyte, and cathode potentials vs. time and space.

   :param sol: A pseudo-2D model solution object.
   :type sol: P2D Solution object

   :returns: *None.*


