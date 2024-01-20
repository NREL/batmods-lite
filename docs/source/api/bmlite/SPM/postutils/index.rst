:py:mod:`bmlite.SPM.postutils`
==============================

.. py:module:: bmlite.SPM.postutils

.. autoapi-nested-parse::

   Post-processing Utilities Module
   --------------------------------
   This module contains all post-processing functions for the SPM package. The
   available post-processing options for a given experiment are specific to that
   experiment. Therefore, not all ``Solution`` classes may have access to all of
   the following functions.



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   bmlite.SPM.postutils.contour
   bmlite.SPM.postutils.contours
   bmlite.SPM.postutils.current
   bmlite.SPM.postutils.debug
   bmlite.SPM.postutils.general
   bmlite.SPM.postutils.post
   bmlite.SPM.postutils.power
   bmlite.SPM.postutils.verify
   bmlite.SPM.postutils.voltage



.. py:function:: contour(ax: object, xlim: list[float], ylim: list[float], z: numpy.ndarray, label: str) -> None


.. py:function:: contours(sol: object) -> None


.. py:function:: current(sol: object) -> None


.. py:function:: debug(sol: object) -> None


.. py:function:: general(sol: object) -> None


.. py:function:: post(sol: object) -> dict


.. py:function:: power(sol: object) -> None


.. py:function:: verify(sol: object) -> None


.. py:function:: voltage(sol: object) -> None


