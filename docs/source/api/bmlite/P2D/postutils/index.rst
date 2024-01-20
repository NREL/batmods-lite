:py:mod:`bmlite.P2D.postutils`
==============================

.. py:module:: bmlite.P2D.postutils

.. autoapi-nested-parse::

   Post-processing Utilities Module
   --------------------------------
   This module contains all post-processing functions for the P2D package. The
   available post-processing options for a given experiment are specific to that
   experiment. Therefore, not all ``Solution`` classes may have access to all of
   the following functions.



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   bmlite.P2D.postutils.contour
   bmlite.P2D.postutils.contours
   bmlite.P2D.postutils.debug
   bmlite.P2D.postutils.general
   bmlite.P2D.postutils.postvars
   bmlite.P2D.postutils.verify



.. py:function:: contour(ax: object, xlim: list[float], ylim: list[float], z: numpy.ndarray, label: str) -> None


.. py:function:: contours(sol: object) -> None


.. py:function:: debug(sol: object) -> None


.. py:function:: general(sol: object) -> None


.. py:function:: postvars(sol: object) -> dict


.. py:function:: verify(sol: object) -> None


