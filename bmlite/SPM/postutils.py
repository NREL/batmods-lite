"""
Post-processing Utilities Module
--------------------------------
This module contains all post-processing functions for the SPM package. The
available post-processing options for a given experiment are specific to that
experiment. Therefore, not all ``Solution`` classes may have access to all of
the following functions.
"""


def post(sol: object) -> dict:
    import numpy as np
    from scipy.integrate import cumtrapz

    from .dae import residuals

    # Pull sim and exp from sol
    sim, exp = sol._sim, sol._exp.copy()

    # Extract desired variables for each time
    res = np.zeros_like(sol.y)

    sdot_an = np.zeros_like(sol.t)
    sdot_ca = np.zeros_like(sol.t)

    i_ext = np.zeros_like(sol.t)

    # Turn on output from residuals
    sim._flags['post'] = True

    for i, t in enumerate(sol.t):
        sv, svdot = sol.y[i, :], sol.ydot[i, :]

        output = residuals(t, sv, svdot, np.zeros_like(sv), (sim, exp))
        res[i, :], sdot_an[i], sdot_ca[i] = output

        i_ext[i] = exp['i_ext']

    # Turn off output from residuals
    sim._flags['post'] = False

    # Areal capacity [A*h/m^2]
    cap_m2 = np.abs(np.hstack([0., cumtrapz(i_ext, sol.t / 3600.)]))

    # Store outputs
    postvars = {}
    postvars['res'] = res

    postvars['sdot_an'] = sdot_an
    postvars['sdot_ca'] = sdot_ca

    postvars['i_ext'] = i_ext
    postvars['A*h/m^2'] = cap_m2

    return postvars


def current(sol: object, ax: object = None) -> None:
    import matplotlib.pyplot as plt

    from ..plotutils import format_ticks, show

    if len(sol.postvars) == 0:
        sol.post()

    if ax is None:
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=[5, 3.5])

    ax.set_xlabel(r'$t$ [s]')
    ax.set_ylabel(r'Current density, $i_{\rm ext}$ [A/m$^2$]')

    ax.plot(sol.t, sol.postvars['i_ext'], '-k')
    format_ticks(ax)

    if ax is None:
        show(fig)


def voltage(sol: object, ax: object = None) -> None:
    import matplotlib.pyplot as plt

    from ..plotutils import format_ticks, show

    if ax is None:
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=[5, 3.5])

    ax.set_xlabel(r'$t$ [s]')
    ax.set_ylabel('Cell voltage [V]')

    sim = sol._sim

    ax.plot(sol.t, sol.y[:, sim.ca.ptr['phi_ed']], '-k')
    format_ticks(ax)

    if ax is None:
        show(fig)


def power(sol: object, ax: object = None) -> None:
    import matplotlib.pyplot as plt

    from ..plotutils import format_ticks, show

    if len(sol.postvars) == 0:
        sol.post()

    if ax is None:
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=[5, 3.5])

    ax.set_xlabel(r'$t$ [s]')
    ax.set_ylabel(r'Power density, $P_{\rm ext}$ [W/m$^2$]')

    sim = sol._sim

    i_ext = sol.postvars['i_ext']
    V_cell = sol.y[:, sim.ca.ptr['phi_ed']]

    ax.plot(sol.t, i_ext * V_cell, '-k')
    format_ticks(ax)

    if ax is None:
        show(fig)


def IVP(sol: object) -> None:
    import matplotlib.pyplot as plt

    from ..plotutils import show

    fig, ax = plt.subplots(nrows=1, ncols=3, figsize=[15, 3.5])

    current(sol, ax[0])
    voltage(sol, ax[1])
    power(sol, ax[2])

    fig.subplots_adjust(wspace=0.3)
    show(fig)


def potentials(sol: object) -> None:
    import matplotlib.pyplot as plt

    from ..plotutils import format_ticks, show

    an, ca, el = sol._sim.an, sol._sim.ca, sol._sim.el

    fig, ax = plt.subplots(nrows=1, ncols=3, figsize=[15, 3.5])

    ax[0].set_ylabel(r'$\phi_{\rm an}$ [V]')
    ax[1].set_ylabel(r'$\phi_{\rm ca}$ [V]')
    ax[2].set_ylabel(r'$\phi_{\rm el}$ [V]')

    ax[0].plot(sol.t, sol.y[:, an.ptr['phi_ed']], '-C3')
    ax[1].plot(sol.t, sol.y[:, ca.ptr['phi_ed']], '-C2')
    ax[2].plot(sol.t, sol.y[:, el.ptr['phi_el']], '-C0')

    for i in range(3):
        ax[i].set_xlabel(r'$t$ [s]')
        format_ticks(ax[i])

    fig.subplots_adjust(wspace=0.3)
    show(fig)


def intercalation(sol: object) -> None:
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.colors as clrs

    from ..plotutils import format_ticks, show

    # Pull sim and exp from sol
    sim = sol._sim

    # Break inputs into separate objects
    an, ca = sim.an, sim.ca

    # Pull time indices and setup colorbar
    t_inds = np.ceil(np.linspace(0, sol.t.size - 1, 11)).astype(int)
    cmap = plt.get_cmap('jet', len(t_inds))

    norm = clrs.Normalize(vmin=sol.t.min(), vmax=sol.t.max())
    sm = plt.cm.ScalarMappable(cmap='jet', norm=norm)

    # Solid-phase Li intercalation fracs -- anode and cathode
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=[8, 3],
                           layout='constrained')

    ax[0].set_ylabel(r'$X_{\rm Li, an}$ [$-$]')
    ax[1].set_ylabel(r'$X_{\rm Li, ca}$ [$-$]')

    ax[0].text(0.1, 0.1, 'Anode particle', transform=ax[0].transAxes)
    ax[1].text(0.1, 0.1, 'Cathode particle', transform=ax[1].transAxes)

    for i, it in enumerate(t_inds):
        Li_an = sol.y[it, an.r_ptr('Li_ed')]
        ax[0].plot(an.r * 1e6, Li_an, color=cmap(i))

        Li_ca = sol.y[it, ca.r_ptr('Li_ed')]
        ax[1].plot(ca.r * 1e6, Li_ca, color=cmap(i))

    cb = plt.colorbar(sm, ax=ax[1], ticks=sol.t[t_inds])
    cb.set_label(r'$t$ [s]')

    ax[0].set_xlim([0., an.R_s * 1e6])
    ax[1].set_xlim([0., ca.R_s * 1e6])

    for i in range(2):
        ax[i].set_xlabel(r'$r$ [$\mu$m]')
        ax[i].set_ylim([0., 1.05])
        format_ticks(ax[i])

    show(fig)


def contours(sol: object) -> None:
    import matplotlib.pyplot as plt

    from ..plotutils import contour, show

    # Get needed domains
    an, ca = sol._sim.an, sol._sim.ca

    # Make figure
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=[6, 3.5])

    # Li concentrations in anode [kmol/m^3]
    xlims = [an.rm[0] * 1e6, an.rp[-1] * 1e6]
    ylims = [sol.t.min(), sol.t.max()]
    z = sol.y[:, an.r_ptr('Li_ed')] * an.Li_max

    contour(ax[0], xlims, ylims, z, r'[kmol/m$^3$]')

    ax[0].set_ylabel(r'$t$ [s]')

    ax[0].set_xlabel(r'$r$ [$\mu$m]')
    ax[0].set_title(r'$C_{\rm s, an}$')

    # Li concentrations in cathode [kmol/m^3]
    xlims = [ca.rm[0] * 1e6, ca.rp[-1] * 1e6]
    ylims = [sol.t.min(), sol.t.max()]
    z = sol.y[:, ca.r_ptr('Li_ed')] * ca.Li_max

    contour(ax[1], xlims, ylims, z, r'[kmol/m$^3$]')

    ax[1].set_xlabel(r'$r$ [$\mu$m]')
    ax[1].set_title(r'$C_{\rm s, ca}$')

    # Adjust spacing
    fig.subplots_adjust(wspace=0.8, hspace=0.2)
    show(fig)
