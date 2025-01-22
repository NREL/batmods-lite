"""
Post-processing Utilities
-------------------------
This module contains all post-processing functions for the P2D package. The
available post-processing options for a given experiment are specific to that
experiment. Therefore, not all ``Solution`` classes may have access to all of
the following functions.

"""

from typing import TypeVar

import numpy as np
import matplotlib.pyplot as plt

from ._solutions import BaseSolution

Solution = TypeVar('Solution', bound='BaseSolution')


def post(soln: Solution) -> dict:
    """
    Run post processing to determine secondary variables.

    Parameters
    ----------
    soln : Solution
        A pseudo-2D model solution object.

    Returns
    -------
    postvars : dict
        Post processed variables, as described below.

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

    See also
    --------
    ~bmlite.P2D.solutions.StepSolution
    ~bmlite.P2D.solutions.CycleSolution

    """

    from .dae import residuals

    # Pull sim and exp from soln
    sim = soln._sim
    step = {
        'mode': 'post',
        'units': 'post',
        'value': 'post',
    }

    # Get needed domains
    an, sep, ca = sim.an, sim.sep, sim.ca

    # Extract desired variables for each time
    div_i_an = np.zeros([soln.t.size, an.Nx])
    div_i_sep = np.zeros([soln.t.size, sep.Nx])
    div_i_ca = np.zeros([soln.t.size, ca.Nx])

    sdot_an = np.zeros([soln.t.size, an.Nx])
    sdot_ca = np.zeros([soln.t.size, ca.Nx])

    sum_ip = np.zeros([soln.t.size, an.Nx + sep.Nx + ca.Nx])
    i_el_x = np.zeros([soln.t.size, an.Nx + sep.Nx + ca.Nx + 1])

    # Turn on output from residuals
    for i, t in enumerate(soln.t):
        sv, svdot = soln.y[i, :], soln.yp[i, :]

        output = residuals(t, sv, svdot, np.zeros_like(sv), (sim, step))

        (div_i_an[i, :], div_i_sep[i, :], div_i_ca[i, :], sdot_an[i, :],
         sdot_ca[i, :], sum_ip[i, :], i_el_x[i, :]) = output

    # Store outputs
    postvars = {
        'div_i_an': div_i_an,
        'div_i_sep': div_i_sep,
        'div_i_ca': div_i_ca,
        'sdot_an': sdot_an,
        'sdot_ca': sdot_ca,
        'sum_ip': sum_ip,
        'i_el_x': i_el_x,
    }

    return postvars


def _liquid_phase_Li(soln: Solution) -> np.ndarray:
    """
    Calculate the liquid-phase lithium vs. time.

    Parameters
    ----------
    soln : SPM Solution object
        A single particle model solution object.

    Returns
    -------
    Li_ed_0 : float
        liquid-phase lithium [kmol/m^2] based on ``el.Li_0``.
    Li_ed_t : 1D array
        Solution's liquid-phase lithium [kmol/m^2] vs. time [s].

    """

    el, an, sep, ca = soln._sim.el, soln._sim.an, soln._sim.sep, soln._sim.ca

    # Initial total liquid-phase lithium [kmol/m^2]
    Li_el_0 = np.sum(np.hstack([an.eps_el * el.Li_0 * (an.xp - an.xm),
                                sep.eps_el * el.Li_0 * (sep.xp - sep.xm),
                                ca.eps_el * el.Li_0 * (ca.xp - ca.xm)]))

    # Total liquid-phase lithium [kmol/m^2] vs. time [s]
    Li_an = soln.vars['an']['ce']
    Li_sep = soln.vars['sep']['ce']
    Li_ca = soln.vars['ca']['ce']

    Li_el_t = np.sum(np.hstack([an.eps_el * Li_an * (an.xp - an.xm),
                                sep.eps_el * Li_sep * (sep.xp - sep.xm),
                                ca.eps_el * Li_ca * (ca.xp - ca.xm)]), axis=1)

    return Li_el_0, Li_el_t


def _solid_phase_Li(soln: Solution) -> np.ndarray:
    """
    Calculate the solid-phase lithium vs. time.

    Parameters
    ----------
    soln : SPM Solution object
        A single particle model solution object.

    Returns
    -------
    Li_ed_0 : float
        Solid-phase lithium [kmol/m^2] based on ``an.x_0`` and ``ca.x_0``.
    Li_ed_t : 1D array
        Solution's solid-phase lithium [kmol/m^2] vs. time [s].

    """

    an, ca = soln._sim.an, soln._sim.ca

    # Initial total solid-phase lithium [kmol/m^2]
    Li_ed_0 = an.x_0*an.Li_max*an.eps_AM*an.thick \
            + ca.x_0*ca.Li_max*ca.eps_AM*ca.thick

    # Anode/cathode lithium [kmol/m^2] vs. time [s]
    V_an = 4.*np.pi*an.R_s**3 / 3.
    V_ca = 4.*np.pi*ca.R_s**3 / 3.

    dr_an = an.rp - an.rm
    dr_ca = ca.rp - ca.rm

    Li_an = np.zeros_like(soln.t)
    Li_ca = np.zeros_like(soln.t)

    for i in range(soln.t.size):

        Li_ed_xr = soln.vars['an']['cs'][i, :, :]
        Li_ed_x = np.sum(4.*np.pi*an.r**2 * Li_ed_xr*dr_an, axis=1) / V_an

        Li_an[i] = np.sum(an.eps_AM*Li_ed_x*(an.xp - an.xm))

        Li_ed_xr = soln.vars['ca']['cs'][i, :, :]
        Li_ed_x = np.sum(4.*np.pi*ca.r**2 * Li_ed_xr*dr_ca, axis=1) / V_ca

        Li_ca[i] = np.sum(ca.eps_AM*Li_ed_x*(ca.xp - ca.xm))

    # Total solid-phase lithium [kmol/m^2] vs. time [s]
    Li_ed_t = Li_an + Li_ca

    return Li_ed_0, Li_ed_t


def potentials(soln: Solution) -> None:
    """
    Plots anode, electrolyte, and cathode potentials vs. time and space.

    Parameters
    ----------
    soln : P2D Solution object
        A pseudo-2D model solution object.

    Returns
    -------
    None.

    """

    import matplotlib.colors as clrs

    from ..plotutils import format_ticks, show

    sep, ca = soln._sim.sep, soln._sim.ca

    # Pull time indices and setup colorbar
    t_inds = np.ceil(np.linspace(0, soln.t.size - 1, 11)).astype(int)

    norm = clrs.Normalize(vmin=soln.t.min(), vmax=soln.t.max())
    sm = plt.cm.ScalarMappable(cmap='Greys', norm=norm)

    # Phase potentials [V] vs. time [s]
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=[8, 3],
                           layout='constrained')

    cmap = plt.get_cmap('Reds', len(t_inds))
    for i, it in enumerate(t_inds):
        if it != t_inds[-4]:
            label = '__nolabel'
        else:
            label = r'$\phi_{\rm an}$'

        ax.plot(soln.vars['an']['x']*1e6, soln.vars['an']['phis'][it, :],
                color=cmap(i), label=label)

    cmap = plt.get_cmap('Blues', len(t_inds))
    for i, it in enumerate(t_inds):
        if it != t_inds[-4]:
            label = '__nolabel'
        else:
            label = r'$\phi_{\rm el}$'

        ax.plot(soln.vars['el']['x']*1e6, soln.vars['el']['phie'][it, :],
                color=cmap(i), label=label)

    cmap = plt.get_cmap('Greens', len(t_inds))
    for i, it in enumerate(t_inds):
        if it != t_inds[-4]:
            label = '__nolabel'
        else:
            label = r'$\phi_{\rm ca}$'

        ax.plot(soln.vars['ca']['x']*1e6, soln.vars['ca']['phis'][it, :],
                color=cmap(i), label=label)

    cb = plt.colorbar(sm, ax=ax, ticks=soln.t[t_inds])
    cb.set_label(r'$t$ [s]')

    ax.set_xlabel(r'$x$ [$\mu$m]')
    ax.set_ylabel(r'Potentials [V]')

    ax.legend(loc='upper left', frameon=False, borderpad=2)

    ax.set_xlim([0., ca.xp[-1]*1e6])

    ylims = ax.get_ylim()
    ax.set_ylim(ylims)

    ax.vlines(sep.xm[0]*1e6, ylims[0], ylims[1], 'k', linestyles='--')
    ax.vlines(sep.xp[-1]*1e6, ylims[0], ylims[1], 'k', linestyles='--')

    format_ticks(ax)
    show(fig)


def electrolyte(soln: Solution) -> None:
    """
    Plots electrolyte Li-ion concentration profiles vs. time.

    Parameters
    ----------
    soln : P2D Solution object
        A pseudo-2D model solution object.

    Returns
    -------
    None.

    """

    import matplotlib.colors as clrs

    from ..plotutils import format_ticks, show

    sep, ca = soln._sim.sep, soln._sim.ca

    # Pull time indices and setup colorbar
    t_inds = np.ceil(np.linspace(0, soln.t.size - 1, 11)).astype(int)
    cmap = plt.get_cmap('jet', len(t_inds))

    norm = clrs.Normalize(vmin=soln.t.min(), vmax=soln.t.max())
    sm = plt.cm.ScalarMappable(cmap='jet', norm=norm)

    # Electrolyte-phase Li-ion concentration [kmol/m^3]
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=[8, 3],
                           layout='constrained')

    for i, it in enumerate(t_inds):
        ax.plot(soln.vars['el']['x']*1e6, soln.vars['el']['ce'][it, :],
                color=cmap(i))

    cb = plt.colorbar(sm, ax=ax, ticks=soln.t[t_inds])
    cb.set_label(r'$t$ [s]')

    ax.set_xlabel(r'$x$ [$\mu$m]')
    ax.set_ylabel(r'$C_{\rm Li^+}$ [kmol/m$^3$]')

    ax.set_xlim([0., ca.xp[-1] * 1e6])

    ylims = ax.get_ylim()
    ax.set_ylim(ylims)

    ax.vlines(sep.xm[0]*1e6, ylims[0], ylims[1], 'k', linestyles='--')
    ax.vlines(sep.xp[-1]*1e6, ylims[0], ylims[1], 'k', linestyles='--')

    format_ticks(ax)
    show(fig)


def intercalation(soln: Solution) -> None:
    """
    Plots anode and cathode particle intercalation profiles vs. time.

    Parameters
    ----------
    soln : P2D Solution object
        A pseudo-2D model solution object.

    Returns
    -------
    None.

    """

    import matplotlib.colors as clrs

    from ..plotutils import format_ticks, show

    # Pull time indices and setup colorbar
    t_inds = np.ceil(np.linspace(0, soln.t.size - 1, 11)).astype(int)
    cmap = plt.get_cmap('jet', len(t_inds))

    norm = clrs.Normalize(vmin=soln.t.min(), vmax=soln.t.max())
    sm = plt.cm.ScalarMappable(cmap='jet', norm=norm)

    # Solid-phase Li intercalation fracs [-]
    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=[8, 6],
                           layout='constrained')

    ax[0, 0].text(0.1, 0.1, r'$x$ = an/sep', transform=ax[0, 0].transAxes)
    ax[0, 1].text(0.1, 0.1, r'$x$ = sep/ca', transform=ax[0, 1].transAxes)
    ax[1, 0].text(0.1, 0.1, r'$x$ = cc/an', transform=ax[1, 0].transAxes)
    ax[1, 1].text(0.1, 0.1, r'$x$ = ca/cc', transform=ax[1, 1].transAxes)

    for i, it in enumerate(t_inds):
        ax[0, 0].plot(soln.vars['an']['r']*1e6, soln.vars['an']['xs'][it, -1],
                      color=cmap(i))

        ax[0, 1].plot(soln.vars['ca']['r']*1e6, soln.vars['ca']['xs'][it, 0],
                      color=cmap(i))

        ax[1, 0].plot(soln.vars['an']['r']*1e6, soln.vars['an']['xs'][it, 0],
                      color=cmap(i))

        ax[1, 1].plot(soln.vars['ca']['r']*1e6, soln.vars['ca']['xs'][it, -1],
                      color=cmap(i))

    cax = ax.ravel().tolist()
    cb = plt.colorbar(sm, ax=cax, ticks=soln.t[t_inds], aspect=50)
    cb.set_label(r'$t$ [s]')

    for i in range(2):
        ax[i, 0].set_ylabel(r'$X_{\rm Li}$ [$-$]')
        ax[i, 1].set_yticklabels([])

    for j in range(2):
        ax[0, j].set_xticklabels([])
        ax[1, j].set_xlabel(r'$r$ [$\mu$m]')

    for i in range(2):
        for j in range(2):
            ax[i, j].set_ylim([0., 1.05])
            format_ticks(ax[i, j])

    show(fig)


def pixels(soln: Solution) -> None:
    """
    Makes pixel plots for most 2D (space/time) variables.

    Parameters
    ----------
    soln : P2D Solution object
        A pseudo-2D model solution object.

    Returns
    -------
    None.

    """

    from ..plotutils import pixel, show

    # Get needed domains
    an, ca = soln._sim.an, soln._sim.ca

    # Make figure
    fig, ax = plt.subplots(nrows=3, ncols=3, figsize=[8, 10],
                           layout='constrained')

    # Li-ion conc. [kmol/m^3]
    xlims = [an.xm[0] * 1e6, ca.xp[-1] * 1e6]
    ylims = [soln.t.min(), soln.t.max()]
    z = soln.vars['el']['ce']

    pixel(ax[0, 0], xlims, ylims, z, r'[kmol/m$^3$]')

    ax[0, 0].set_ylabel(r'$t$ [s]')
    ax[0, 0].set_title(r'$C_{\rm Li+}$')

    # Electrolyte potential [V]
    xlims = [an.xm[0] * 1e6, ca.xp[-1] * 1e6]
    ylims = [soln.t.min(), soln.t.max()]
    z = soln.vars['el']['phie']

    pixel(ax[0, 1], xlims, ylims, z, r'[V]')

    ax[0, 1].set_yticks([])
    ax[0, 1].set_title(r'$\phi_{\rm el}$')

    # Ionic current [A/m^2]
    xlims = [an.xm[0] * 1e6, ca.xp[-1] * 1e6]
    ylims = [soln.t.min(), soln.t.max()]
    z = soln.vars['el']['ie']

    pixel(ax[0, 2], xlims, ylims, z, r'[A/m$^2$]')

    ax[0, 2].set_yticks([])
    ax[0, 2].set_title(r'$i_{\rm el}$')

    # Surface conc. for anode [kmol/m^3]
    xlims = [an.x[0] * 1e6, an.x[-1] * 1e6]
    ylims = [soln.t.min(), soln.t.max()]
    z = soln.vars['an']['cs'][:, :, -1]

    pixel(ax[1, 0], xlims, ylims, z, r'[kmol/m$^3$]')

    ax[1, 0].set_ylabel(r'$t$ [s]')
    ax[1, 0].set_title(r'$C_{\rm s, an}$')

    # Anode potential [mV]
    xlims = [an.x[0] * 1e6, an.x[-1] * 1e6]
    ylims = [soln.t.min(), soln.t.max()]
    z = soln.vars['an']['phis']

    pixel(ax[1, 1], xlims, ylims, z * 1e3, r'[mV]')

    ax[1, 1].set_yticks([])
    ax[1, 1].set_title(r'$\phi_{\rm s, an}$')

    # Faradaic current in anode [kmol/m^2/s]
    xlims = [an.x[0] * 1e6, an.x[-1] * 1e6]
    ylims = [soln.t.min(), soln.t.max()]
    z = soln.vars['an']['sdot']

    pixel(ax[1, 2], xlims, ylims, z, r'[kmol/m$^2$/s]')

    ax[1, 2].set_yticks([])
    ax[1, 2].set_title(r'$j_{\rm Far, an}$')

    # Surface conc. for cathode [kmol/m^3]
    xlims = [ca.x[0] * 1e6, ca.x[-1] * 1e6]
    ylims = [soln.t.min(), soln.t.max()]
    z = soln.vars['ca']['cs'][:, :, -1]

    pixel(ax[2, 0], xlims, ylims, z, r'[kmol/m$^3$]')

    ax[2, 0].set_ylabel(r'$t$ [s]')
    ax[2, 0].set_xlabel(r'$x$ [$\mu$m]')
    ax[2, 0].set_title(r'$C_{\rm s, ca}$')

    # Cathode potential [V]
    xlims = [ca.x[0] * 1e6, ca.x[-1] * 1e6]
    ylims = [soln.t.min(), soln.t.max()]
    z = soln.vars['ca']['phis']

    pixel(ax[2, 1], xlims, ylims, z, r'[V]')

    ax[2, 1].set_yticks([])
    ax[2, 1].set_xlabel(r'$x$ [$\mu$m]')
    ax[2, 1].set_title(r'$\phi_{\rm s, ca}$')

    # Faradaic current in cathode [kmol/m^2/s]
    xlims = [ca.x[0] * 1e6, ca.x[-1] * 1e6]
    ylims = [soln.t.min(), soln.t.max()]
    z = soln.vars['ca']['sdot']

    pixel(ax[2, 2], xlims, ylims, z, r'[kmol/m$^2$/s]')

    ax[2, 2].set_yticks([])
    ax[2, 2].set_xlabel(r'$x$ [$\mu$m]')
    ax[2, 2].set_title(r'$j_{\rm Far, ca}$')

    # Adjust spacing
    fig.get_layout_engine().set(hspace=0.1, wspace=0.1)
    show(fig)
