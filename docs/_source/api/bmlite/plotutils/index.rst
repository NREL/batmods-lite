:py:mod:`bmlite.plotutils`
==========================

.. py:module:: bmlite.plotutils

.. autoapi-nested-parse::

   .. rubric:: Plotting Utilities

   A module with functions for plotting data and formatting figures. Functions
   here are generally useful for all models in BatMods-lite. More specific plots
   are written within the ``postutils`` modules of their respective model.



Package Contents
----------------


Functions
~~~~~~~~~

.. autoapisummary::

   bmlite.plotutils.format_ticks
   bmlite.plotutils.pixel
   bmlite.plotutils.show



.. py:function:: format_ticks(ax)

   Formats an ``axis`` object by adjusting the ticks.

   Specifically, the top and right ticks are added, minor ticks are turned
   on, and all ticks are set to face inward.

   :param ax: An ``axis`` instance from a ``matplotlib`` figure.
   :type ax: object

   :returns: *None.*


.. py:function:: pixel(ax, xlims, ylims, z, cblabel)

   Fill an axis instance with a pixel plot defined by the inputs.

   :param ax: An ``axis`` instance from a ``matplotlib`` figure.
   :type ax: object
   :param xlims: Limits for the x-axis [x_low, x_high].
   :type xlims: list[float]
   :param ylims: Limits for the y-axis [y_low, y_high].
   :type ylims: list[float]
   :param z: Data to plot against x and y. ``z[0, 0]`` corresponds to x_low, y_low,
             and ``z[-1, -1]`` corresponds to x_high, y_high.
   :type z: 2D array
   :param cblabel: The colorbar label.
   :type cblabel: str

   :returns: *None.*


.. py:function:: show(fig)

   Display a figure according to the backend.

   :param fig: A ``fig`` instance from a ``matplotlib`` figure.
   :type fig: object

   :returns: *None.*


