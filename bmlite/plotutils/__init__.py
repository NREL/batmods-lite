"""
Plotting Utilities
------------------
A module with functions for plotting data and formatting figures. Functions
here are generally useful for all models in BATMODS-lite. More specific plots
are written within the ``postutils`` modules of their respective model.
"""

from numpy import ndarray as _ndarray


def show(fig: object) -> None:
    """
    Display a figure according to the backend.

    Parameters
    ----------
    fig : object
        A ``fig`` instance from a ``matplotlib`` figure.

    Returns
    -------
    None.
    """

    from matplotlib import get_backend

    if 'inline' not in get_backend():
        fig.show()


def format_ticks(ax: object) -> None:
    """
    Formats an ``axis`` object by adjusting the ticks.

    Specifically, the top and right ticks are added, minor ticks are turned
    on, and all ticks are set to face inward.

    Parameters
    ----------
    ax : object
        An ``axis`` instance from a ``matplotlib`` figure.

    Returns
    -------
    None.
    """

    from matplotlib.ticker import AutoMinorLocator

    if ax.get_xaxis().get_scale() != 'log':
        ax.xaxis.set_minor_locator(AutoMinorLocator())

    if ax.get_yaxis().get_scale() != 'log':
        ax.yaxis.set_minor_locator(AutoMinorLocator())

    ax.tick_params(axis='x', top=True, which='both', direction='in')
    ax.tick_params(axis='y', right=True, which='both', direction='in')


def pixel(ax: object, xlims: list[float], ylims: list[float], z: _ndarray,
            cblabel: str) -> None:
    """
    Fill an axis instance with a pixel plot defined by the inputs.

    Parameters
    ----------
    ax : object
        An ``axis`` instance from a ``matplotlib`` figure.

    xlims : list[float]
        Limits for the x-axis [x_low, x_high].

    ylims : list[float]
        Limits for the y-axis [y_low, y_high].

    z : 2D array
        Data to plot against x and y. ``z[0, 0]`` corresponds to x_low, y_low,
        and ``z[-1, -1]`` corresponds to x_high, y_high.

    cblabel : str
        The colorbar label.

    Returns
    -------
    None.
    """

    import matplotlib.pyplot as plt

    ax.set_xticks([])

    cmap = plt.cm.viridis

    im = ax.imshow(z, cmap=cmap, aspect='auto', vmin=z.min(), vmax=z.max(),
                   extent=[xlims[0], xlims[1], ylims[1], ylims[0]],
                   interpolation='nearest')

    cb = plt.colorbar(im, ax=ax)
    cb.ax.yaxis.set_offset_position('left')
    cb.set_label(cblabel)
