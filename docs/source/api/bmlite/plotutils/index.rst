:py:mod:`bmlite.plotutils`
==========================

.. py:module:: bmlite.plotutils


Package Contents
----------------


Functions
~~~~~~~~~

.. autoapisummary::

   bmlite.plotutils.format_lims
   bmlite.plotutils.format_ticks



.. py:function:: format_lims(ax: object) -> None

   Formats an ``axis`` object by adjusting the limits.

   Specifically, sets the x and y limits such that there is 10% white space
   all the way around the perimeter. In the case that x or y is constant, the
   default matplotlib behavior is adopted for the limits on that respective
   variable.

   :param ax: An ``axis`` instance from a ``matplotlib`` figure.
   :type ax: object

   :returns: *None.*


.. py:function:: format_ticks(ax: object) -> None

   Formats an ``axis`` object by adjusting the ticks.

   Specifically, the top and right ticks are added, minor ticks are turned on,
   and all ticks are set to face inward.

   :param ax: An ``axis`` instance from a ``matplotlib`` figure.
   :type ax: object

   :returns: *None.*


