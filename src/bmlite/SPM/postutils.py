"""
Post-processing Utilities
-------------------------
This module contains all post-processing functions for the SPM package. The
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
        A single particle model solution object.

    Returns
    -------
    postvars : dict
        Post processed variables, as described below.

        ========= ========================================================
        Key       Value [units] (*type*)
        ========= ========================================================
        sdot_an   anode Faradaic current at t [kmol/m^2/s] (*1D array*)
        sdot_ca   cathode Faradaic current at t [kmol/m^2/s] (*1D array*)
        ========= ========================================================

    See also
    --------
    ~bmlite.SPM.solutions.StepSolution
    ~bmlite.SPM.solutions.CycleSolution

    """

    from .dae import residuals

    # Pull sim and exp from sol
    sim = soln._sim
    step = {
        'mode': 'post',
        'units': 'post',
        'value': 'post',
    }

    # Extract desired variables for each time
    sdot_an = np.zeros_like(soln.t)
    sdot_ca = np.zeros_like(soln.t)

    # Turn on output from residuals
    for i, t in enumerate(soln.t):
        sv, svdot = soln.y[i, :], soln.yp[i, :]

        output = residuals(t, sv, svdot, np.zeros_like(sv), (sim, step))
        sdot_an[i], sdot_ca[i] = output

    # Store outputs
    postvars = {
        'sdot_an': sdot_an,
        'sdot_ca': sdot_ca,
    }

    return postvars


def _solid_phase_Li(soln: Solution) -> np.array:
    """
    Calculate the solid-phase lithium vs. time.

    Parameters
    ----------
    soln : Solution
        A single particle model solution object.

    Returns
    -------
    Li_ed_0 : float
        Solid-phase lithium [kmol/m^2] based on ``an.x_0`` and ``ca.x_0``.
    Li_ed_t : 1D array
        Solution's solid-phase lithium [kmol/m^2] vs. time [s].

    See also
    --------
    ~bmlite.SPM.solutions.StepSolution
    ~bmlite.SPM.solutions.CycleSolution

    """

    from ..mathutils import int_r

    an, ca = soln._sim.an, soln._sim.ca

    # Initial total solid-phase lithium [kmol/m^2]
    Li_ed_0 = an.x_0*an.Li_max*an.eps_AM*an.thick \
            + ca.x_0*ca.Li_max*ca.eps_AM*ca.thick

    # Anode/cathode lithium [kmol/m^2] vs. time [s]
    V_an = 4.*np.pi*an.R_s**3 / 3.
    V_ca = 4.*np.pi*ca.R_s**3 / 3.

    Li_an = np.zeros_like(soln.t)
    Li_ca = np.zeros_like(soln.t)

    for i in range(soln.t.size):
        Li_an[i] = an.thick*an.eps_AM / V_an \
                 * int_r(an.rm, an.rp, soln.vars['an']['cs'][i, :])
        Li_ca[i] = ca.thick*ca.eps_AM / V_ca \
                 * int_r(ca.rm, ca.rp, soln.vars['ca']['cs'][i, :])

    # Total solid-phase lithium [kmol/m^2] vs. time [s]
    Li_ed_t = Li_an + Li_ca

    return Li_ed_0, Li_ed_t


def potentials(soln: Solution) -> None:
    """
    Plots anode, electrolyte, and cathode potentials vs. time.

    Parameters
    ----------
    soln : Solution
        A single particle model solution object.

    Returns
    -------
    None.

    See also
    --------
    ~bmlite.SPM.solutions.StepSolution
    ~bmlite.SPM.solutions.CycleSolution

    """

    from ..plotutils import format_ticks, show

    fig, ax = plt.subplots(nrows=1, ncols=3, figsize=[12, 3],
                           layout='constrained')

    ax[0].set_ylabel(r'$\phi_{\rm an}$ [V]')
    ax[1].set_ylabel(r'$\phi_{\rm ca}$ [V]')
    ax[2].set_ylabel(r'$\phi_{\rm el}$ [V]')

    ax[0].plot(soln.vars['time_s'], soln.vars['an']['phis'], '-C3')
    ax[1].plot(soln.vars['time_s'], soln.vars['ca']['phis'], '-C2')
    ax[2].plot(soln.vars['time_s'], soln.vars['el']['phie'], '-C0')

    for i in range(3):
        ax[i].set_xlabel(r'$t$ [s]')
        format_ticks(ax[i])

    fig.get_layout_engine().set(wspace=0.1)
    show(fig)


def intercalation(soln: Solution) -> None:
    """
    Plots anode and cathode particle intercalation profiles vs. time.

    Parameters
    ----------
    soln : Solution
        A single particle model solution object.

    Returns
    -------
    None.

    See also
    --------
    ~bmlite.SPM.solutions.StepSolution
    ~bmlite.SPM.solutions.CycleSolution

    """

    import matplotlib.colors as clrs

    from ..plotutils import format_ticks, show

    # Pull sim and exp from sol
    sim = soln._sim

    # Break inputs into separate objects
    an, ca = sim.an, sim.ca

    # Pull time indices and setup colorbar
    t_inds = np.ceil(np.linspace(0, soln.t.size - 1, 11)).astype(int)
    cmap = plt.get_cmap('jet', len(t_inds))

    norm = clrs.Normalize(vmin=soln.t.min(), vmax=soln.t.max())
    sm = plt.cm.ScalarMappable(cmap='jet', norm=norm)

    # Solid-phase Li intercalation fracs -- anode and cathode
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=[8, 3],
                           layout='constrained')

    ax[0].set_ylabel(r'$X_{\rm Li}$ [$-$]')
    ax[1].set_yticklabels([])

    ax[0].text(0.1, 0.1, 'Anode particle', transform=ax[0].transAxes)
    ax[1].text(0.1, 0.1, 'Cathode particle', transform=ax[1].transAxes)

    for i, it in enumerate(t_inds):
        Li_an = soln.vars['an']['xs'][it, :]
        ax[0].plot(an.r*1e6, Li_an, color=cmap(i))

        Li_ca = soln.vars['ca']['xs'][it, :]
        ax[1].plot(ca.r*1e6, Li_ca, color=cmap(i))

    cb = plt.colorbar(sm, ax=ax[1], ticks=soln.t[t_inds])
    cb.set_label(r'$t$ [s]')

    ax[0].set_xlim([0., an.R_s*1e6])
    ax[1].set_xlim([0., ca.R_s*1e6])

    for i in range(2):
        ax[i].set_xlabel(r'$r$ [$\mu$m]')
        ax[i].set_ylim([0., 1.05])
        format_ticks(ax[i])

    show(fig)


def pixels(soln: Solution) -> None:
    """
    Makes pixel plots for most 2D (space/time) variables.

    Parameters
    ----------
    soln : SPM Solution object
        A single particle model solution object.

    Returns
    -------
    None.

    See also
    --------
    ~bmlite.SPM.solutions.StepSolution
    ~bmlite.SPM.solutions.CycleSolution

    """

    from ..plotutils import pixel, show

    # Get needed domains
    an, ca = soln._sim.an, soln._sim.ca

    # Make figure
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=[5.5, 3.25],
                           layout='constrained')

    # Li concentrations in anode [kmol/m^3]
    xlims = [an.rm[0]*1e6, an.rp[-1]*1e6]
    ylims = [soln.t.min(), soln.t.max()]
    z = soln.vars['an']['cs']

    pixel(ax[0], xlims, ylims, z, r'[kmol/m$^3$]')

    ax[0].set_ylabel(r'$t$ [s]')

    ax[0].set_xlabel(r'$r$ [$\mu$m]')
    ax[0].set_title(r'$C_{\rm s, an}$')

    # Li concentrations in cathode [kmol/m^3]
    xlims = [ca.rm[0]*1e6, ca.rp[-1]*1e6]
    ylims = [soln.t.min(), soln.t.max()]
    z = soln.vars['ca']['cs']

    pixel(ax[1], xlims, ylims, z, r'[kmol/m$^3$]')

    ax[1].set_yticks([])
    ax[1].set_xlabel(r'$r$ [$\mu$m]')
    ax[1].set_title(r'$C_{\rm s, ca}$')

    # Adjust spacing
    fig.get_layout_engine().set(hspace=0.1, wspace=0.1)
    show(fig)
