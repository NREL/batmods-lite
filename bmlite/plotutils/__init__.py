"""
Plotting Utilities
------------------
A module with functions for plotting data and formatting figures. Functions
here are generally usefull for all models in BatMods-lite. More specific plots
are written within the ``postutils`` modules of their respective model.
"""

from numpy import ndarray as _ndarray


def show(fig: object) -> None:
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


def contour(ax: object, xlim: list[float], ylim: list[float], z: _ndarray,
            label: str) -> None:

    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    ax.set_xticks([])

    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='10%', pad=0.2)
    cax.tick_params(direction='in')

    cmap = plt.cm.viridis

    im = ax.imshow(z, cmap=cmap, aspect='auto', vmin=z.min(), vmax=z.max(),
                   extent=[xlim[0], xlim[1], ylim[1], ylim[0]],
                   interpolation='nearest')

    fig = ax.get_figure()

    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label(label)
