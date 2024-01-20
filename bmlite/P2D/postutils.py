"""
Post-processing Utilities Module
--------------------------------
This module contains all post-processing functions for the P2D package. The
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

    # Get needed domains
    an, sep, ca = sim.an, sim.sep, sim.ca

    # Extract desired variables for each time
    res = np.zeros_like(sol.y)

    sdot_an = np.zeros([sol.t.size, an.Nx])
    sdot_ca = np.zeros([sol.t.size, ca.Nx])

    sum_ip = np.zeros([sol.t.size, an.Nx + sep.Nx + ca.Nx])
    i_el_x = np.zeros([sol.t.size, an.Nx + sep.Nx + ca.Nx + 1])

    i_ext = np.zeros_like(sol.t)

    # Turn on output from residuals
    sim._flags['post'] = True

    for i, t in enumerate(sol.t):
        sv, svdot = sol.y[i, :], sol.ydot[i, :]

        (res[i, :], sdot_an[i, :], sdot_ca[i, :], sum_ip[i, :],
         i_el_x[i, :]) = residuals(t, sv, svdot, np.zeros_like(sv), (sim, exp))

        i_ext[i] = exp['i_ext']

    div_i_an = res[:, an.x_ptr('phi_ed')] + res[:, an.x_ptr('phi_el')]
    div_i_sep = res[:, sep.x_ptr('phi_el')]
    div_i_ca = res[:, ca.x_ptr('phi_ed')] + res[:, ca.x_ptr('phi_el')]

    # Turn off output from residuals
    sim._flags['post'] = False

    # Areal capacity [A*h/m^2]
    cap_m2 = np.abs(np.hstack([0., cumtrapz(i_ext, sol.t / 3600.)]))

    # Store outputs
    postvars = {}
    postvars['res'] = res

    postvars['div_i_an'] = div_i_an
    postvars['div_i_sep'] = div_i_sep
    postvars['div_i_ca'] = div_i_ca

    postvars['sdot_an'] = sdot_an
    postvars['sdot_ca'] = sdot_ca

    postvars['sum_ip'] = sum_ip
    postvars['i_el_x'] = i_el_x

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

    ax.plot(sol.t, sol.y[:, sim.ca.x_ptr('phi_ed')[-1]], '-k')
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
    V_cell = sol.y[:, sim.ca.x_ptr('phi_ed')[-1]]

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
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.colors as clrs

    from ..plotutils import format_ticks, show

    an, sep, ca = sol._sim.an, sol._sim.sep, sol._sim.ca

    # Pull time indices and setup colorbar
    t_inds = np.ceil(np.linspace(0, sol.t.size - 1, 11)).astype(int)

    norm = clrs.Normalize(vmin=sol.t.min(), vmax=sol.t.max())
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

        ax.plot(an.x * 1e6, sol.y[it, an.x_ptr('phi_ed')], color=cmap(i),
                label=label)

    x_el = np.hstack([an.x, sep.x, ca.x])
    phi_el = np.hstack([sol.y[:, an.x_ptr('phi_el')],
                        sol.y[:, sep.x_ptr('phi_el')],
                        sol.y[:, ca.x_ptr('phi_el')]])

    cmap = plt.get_cmap('Blues', len(t_inds))
    for i, it in enumerate(t_inds):
        if it != t_inds[-4]:
            label = '__nolabel'
        else:
            label = r'$\phi_{\rm el}$'

        ax.plot(x_el * 1e6, phi_el[it, :], color=cmap(i), label=label)

    cmap = plt.get_cmap('Greens', len(t_inds))
    for i, it in enumerate(t_inds):
        if it != t_inds[-4]:
            label = '__nolabel'
        else:
            label = r'$\phi_{\rm ca}$'

        ax.plot(ca.x * 1e6, sol.y[it, ca.x_ptr('phi_ed')], color=cmap(i),
                label=label)

    cb = plt.colorbar(sm, ax=ax, ticks=sol.t[t_inds])
    cb.set_label(r'$t$ [s]')

    ax.set_xlabel(r'$x$ [$\mu$m]')
    ax.set_ylabel(r'Potentials [V]')

    ax.legend(loc='upper left', frameon=False, borderpad=2)

    ax.set_xlim([0., ca.xp[-1] * 1e6])

    ylims = ax.get_ylim()
    ax.set_ylim(ylims)

    ax.vlines(sep.xm[0] * 1e6, ylims[0], ylims[1], 'k', linestyles='--')
    ax.vlines(sep.xp[-1] * 1e6, ylims[0], ylims[1], 'k', linestyles='--')

    format_ticks(ax)
    show(fig)


def electrolyte(sol: object) -> None:
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.colors as clrs

    from ..plotutils import format_ticks, show

    an, sep, ca = sol._sim.an, sol._sim.sep, sol._sim.ca

    # Pull time indices and setup colorbar
    t_inds = np.ceil(np.linspace(0, sol.t.size - 1, 11)).astype(int)
    cmap = plt.get_cmap('jet', len(t_inds))

    norm = clrs.Normalize(vmin=sol.t.min(), vmax=sol.t.max())
    sm = plt.cm.ScalarMappable(cmap='jet', norm=norm)

    # Electrolyte-phase Li-ion concentration [kmol/m^3]
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=[8, 3],
                           layout='constrained')

    x_el = np.hstack([an.x, sep.x, ca.x])

    Li_el = np.hstack([sol.y[:, an.x_ptr('Li_el')],
                       sol.y[:, sep.x_ptr('Li_el')],
                       sol.y[:, ca.x_ptr('Li_el')]])

    for i, it in enumerate(t_inds):
        ax.plot(x_el * 1e6, Li_el[it, :], color=cmap(i))

    cb = plt.colorbar(sm, ax=ax, ticks=sol.t[t_inds])
    cb.set_label(r'$t$ [s]')

    ax.set_xlabel(r'$x$ [$\mu$m]')
    ax.set_ylabel(r'$C_{\rm Li^+}$ [kmol/m$^3$]')

    ax.set_xlim([0., ca.xp[-1] * 1e6])

    ylims = ax.get_ylim()
    ax.set_ylim(ylims)

    ax.vlines(sep.xm[0] * 1e6, ylims[0], ylims[1], 'k', linestyles='--')
    ax.vlines(sep.xp[-1] * 1e6, ylims[0], ylims[1], 'k', linestyles='--')

    format_ticks(ax)
    show(fig)


def intercalation(sol: object) -> None:
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.colors as clrs

    from ..plotutils import format_ticks, show

    an, ca = sol._sim.an, sol._sim.ca

    # Pull time indices and setup colorbar
    t_inds = np.ceil(np.linspace(0, sol.t.size - 1, 11)).astype(int)
    cmap = plt.get_cmap('jet', len(t_inds))

    norm = clrs.Normalize(vmin=sol.t.min(), vmax=sol.t.max())
    sm = plt.cm.ScalarMappable(cmap='jet', norm=norm)

    # Solid-phase Li intercalation fracs [-]
    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=[8, 6])

    ax[0, 0].text(0.1, 0.1, r'$x$ = an/sep', transform=ax[0, 0].transAxes)
    ax[0, 1].text(0.1, 0.1, r'$x$ = sep/ca', transform=ax[0, 1].transAxes)
    ax[1, 0].text(0.1, 0.1, r'$x$ = cc/an', transform=ax[1, 0].transAxes)
    ax[1, 1].text(0.1, 0.1, r'$x$ = ca/cc', transform=ax[1, 1].transAxes)

    for i, it in enumerate(t_inds):
        an_sep = an.x_ptr('Li_ed')[-1]
        ax[0, 0].plot(an.r * 1e6, sol.y[it, an_sep:an_sep + an.Nr],
                      color=cmap(i))

        sep_ca = ca.x_ptr('Li_ed')[0]
        ax[0, 1].plot(ca.r * 1e6, sol.y[it, sep_ca:sep_ca + ca.Nr],
                      color=cmap(i))

        cc_an = an.x_ptr('Li_ed')[0]
        ax[1, 0].plot(an.r * 1e6, sol.y[it, cc_an:cc_an + an.Nr],
                      color=cmap(i))

        ca_cc = ca.x_ptr('Li_ed')[-1]
        ax[1, 1].plot(ca.r * 1e6, sol.y[it, ca_cc:ca_cc + ca.Nr],
                      color=cmap(i))

    cax = ax.ravel().tolist()
    cb = plt.colorbar(sm, ax=cax, ticks=sol.t[t_inds], aspect=50)
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


def contours(sol: object) -> None:
    import numpy as np
    import matplotlib.pyplot as plt

    from ..plotutils import contour, show

    # Check for postvars
    if len(sol.postvars) == 0:
        sol.postvars = post(sol)

    postvars = sol.postvars

    # Get needed domains
    an, sep, ca = sol._sim.an, sol._sim.sep, sol._sim.ca

    # Make figure
    fig, ax = plt.subplots(nrows=3, ncols=3, figsize=[8, 10])

    # Li-ion conc. [kmol/m^3]
    xlims = [an.xm[0] * 1e6, ca.xp[-1] * 1e6]
    ylims = [sol.t.min(), sol.t.max()]
    z = np.hstack([sol.y[:, an.x_ptr('Li_el')],
                   sol.y[:, sep.x_ptr('Li_el')],
                   sol.y[:, ca.x_ptr('Li_el')]])

    contour(ax[0, 0], xlims, ylims, z, r'[kmol/m$^3$]')

    ax[0, 0].set_ylabel(r'$t$ [s]')
    ax[0, 0].set_title(r'$C_{\rm Li+}$')

    # Electrolyte potential [V]
    xlims = [an.xm[0] * 1e6, ca.xp[-1] * 1e6]
    ylims = [sol.t.min(), sol.t.max()]
    z = np.hstack([sol.y[:, an.x_ptr('phi_el')],
                   sol.y[:, sep.x_ptr('phi_el')],
                   sol.y[:, ca.x_ptr('phi_el')]])

    contour(ax[0, 1], xlims, ylims, z, r'[V]')

    ax[0, 1].set_yticks([])
    ax[0, 1].set_title(r'$\phi_{\rm el}$')

    # Ionic current [A/m^2]
    xlims = [an.xm[0] * 1e6, ca.xp[-1] * 1e6]
    ylims = [sol.t.min(), sol.t.max()]
    z = postvars['i_el_x']

    contour(ax[0, 2], xlims, ylims, z, r'[A/m$^2$]')

    ax[0, 2].set_yticks([])
    ax[0, 2].set_title(r'$i_{\rm el}$')

    # Surface conc. for anode [kmol/m^3]
    xlims = [an.x[0] * 1e6, an.x[-1] * 1e6]
    ylims = [sol.t.min(), sol.t.max()]
    z = sol.y[:, an.x_ptr('Li_ed', an.Nr - 1)] * an.Li_max

    contour(ax[1, 0], xlims, ylims, z, r'[kmol/m$^3$]')

    ax[1, 0].set_ylabel(r'$t$ [s]')
    ax[1, 0].set_title(r'$C_{\rm s, an}$')

    # Anode potential [mV]
    xlims = [an.x[0] * 1e6, an.x[-1] * 1e6]
    ylims = [sol.t.min(), sol.t.max()]
    z = sol.y[:, an.x_ptr('phi_ed')]

    contour(ax[1, 1], xlims, ylims, z * 1e3, r'[mV]')

    ax[1, 1].set_yticks([])
    ax[1, 1].set_title(r'$\phi_{\rm s, an}$')

    # Faradaic current in anode [kmol/m^2/s]
    xlims = [an.x[0] * 1e6, an.x[-1] * 1e6]
    ylims = [sol.t.min(), sol.t.max()]
    z = postvars['sdot_an']

    contour(ax[1, 2], xlims, ylims, z, r'[kmol/m$^2$/s]')

    ax[1, 2].set_yticks([])
    ax[1, 2].set_title(r'$j_{\rm Far, an}$')

    # Surface conc. for cathode [kmol/m^3]
    xlims = [ca.x[0] * 1e6, ca.x[-1] * 1e6]
    ylims = [sol.t.min(), sol.t.max()]
    z = sol.y[:, ca.x_ptr('Li_ed', ca.Nr - 1)] * ca.Li_max

    contour(ax[2, 0], xlims, ylims, z, r'[kmol/m$^3$]')

    ax[2, 0].set_ylabel(r'$t$ [s]')
    ax[2, 0].set_xlabel(r'$x$ [$\mu$m]')
    ax[2, 0].set_title(r'$C_{\rm s, ca}$')

    # Cathode potential [V]
    xlims = [ca.x[0] * 1e6, ca.x[-1] * 1e6]
    ylims = [sol.t.min(), sol.t.max()]
    z = sol.y[:, ca.x_ptr('phi_ed')]

    contour(ax[2, 1], xlims, ylims, z, r'[V]')

    ax[2, 1].set_yticks([])
    ax[2, 1].set_xlabel(r'$x$ [$\mu$m]')
    ax[2, 1].set_title(r'$\phi_{\rm s, ca}$')

    # Faradaic current in cathode [kmol/m^2/s]
    xlims = [ca.x[0] * 1e6, ca.x[-1] * 1e6]
    ylims = [sol.t.min(), sol.t.max()]
    z = postvars['sdot_ca']

    contour(ax[2, 2], xlims, ylims, z, r'[kmol/m$^2$/s]')

    ax[2, 2].set_yticks([])
    ax[2, 2].set_xlabel(r'$x$ [$\mu$m]')
    ax[2, 2].set_title(r'$j_{\rm Far, ca}$')

    # Adjust spacing
    fig.subplots_adjust(wspace=0.7, hspace=0.2)

    show(fig)
