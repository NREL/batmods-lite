"""
DAE Module
----------
This module includes the system of differential algebraic equations (DAE) for
the P2D model. In addition, the ``bandwidth`` function is included in this
module, which helps determine the lower and upper bandwidths of ``residuals``
so the ``'band'`` linear solver option can be used in the ``IDASolver`` class.

"""

import numpy as np

if not hasattr(np, 'concat'):  # pragma: no cover
    np.concat = np.concatenate


def sign(x):
    """Return +1 for x >= 0, -1 for x < 0."""
    return (x >= 0).astype(float) * 2 - 1


def residuals(t: float, sv: np.ndarray, svdot: np.ndarray, res: np.ndarray,
              inputs: tuple[object, dict]) -> None | tuple[np.ndarray]:
    """
    The DAE residuals ``res = M*y' - f(t, y)`` for the P2D model.

    Parameters
    ----------
    t : float
        Value of time [s].
    sv : 1D array
        Solution/state variables at time ``t``.
    svdot : 1D array
        Solution/state variable time derivatives at time ``t``.
    res : 1D array
        An array the same size as ``sv`` and ``svdot``. The values are filled
        in with ``res = M*y' - f(t, y)`` inside this function.
    inputs : (sim : P2D Simulation object, exp : experiment dict)
        The simulation object and experimental details dictionary inputs that
        describe the specific battery and experiment to simulate.

    Returns
    -------
    outputs : tuple[np.ndarray]
        If the experimental step `mode` is set to `post`, then the following
        post-processed variables will be returned in a tuple. Otherwise,
        returns None.

        ========== ======================================================
        Variable   Description [units] (*type*)
        ========== ======================================================
        div_i_an   divergence in anode volume [A/m3] (*1D array*)
        div_i_sep  divergence in separator volume [A/m3] (*1D array*)
        div_i_ca   divergence in cathode volume [A/m3] (*1D array*)
        sdot_an    Li+ production at each x_a [kmol/m3/s] (*1D array*)
        sdot_ca    Li+ production at each x_c [kmol/m3/s] (*1D array*)
        sum_ip     i_ed + i_el at "plus" interfaces [A/m2] (*1D array*)
        i_el_x     i_el at each x interface [A/m2] (*1D array*)
        ========== ======================================================

    """

    from ..mathutils import grad_x, grad_r, div_r

    # Break inputs into separate objects
    sim, exp = inputs

    c, bat, el, an, sep, ca = sim.c, sim.bat, sim.el, sim.an, sim.sep, sim.ca

    # Simulation temperature
    T = bat.temp

    # Organize values from sv
    phi_an = sv[an.x_ptr['phi_ed']]
    Li_el_an = sv[an.x_ptr['Li_el']]
    phi_el_an = sv[an.x_ptr['phi_el']]

    xs_an = sv[an.xr_ptr['Li_ed'].flatten()]
    xs_an = xs_an.reshape([an.Nx, an.Nr])
    Li_an = xs_an*an.Li_max

    if 'Hysteresis' in an._submodels:
        hyst_an = sv[an.x_ptr['hyst']]
        Hyst_an = an.get_Mhyst(xs_an[:, -1])*hyst_an
    else:
        Hyst_an = 0.

    Li_el_sep = sv[sep.x_ptr['Li_el']]
    phi_el_sep = sv[sep.x_ptr['phi_el']]

    phi_ca = sv[ca.x_ptr['phi_ed']]
    Li_el_ca = sv[ca.x_ptr['Li_el']]
    phi_el_ca = sv[ca.x_ptr['phi_el']]

    xs_ca = sv[ca.xr_ptr['Li_ed'].flatten()]
    xs_ca = xs_ca.reshape([ca.Nx, ca.Nr])
    Li_ca = xs_ca*ca.Li_max

    if 'Hysteresis' in ca._submodels:
        hyst_ca = sv[ca.x_ptr['hyst']]
        Hyst_ca = ca.get_Mhyst(xs_ca[:, -1])*hyst_ca
    else:
        Hyst_ca = 0.

    # Pre-calculate ln(Li_el)
    ln_Li_an = np.log(Li_el_an)
    ln_Li_sep = np.log(Li_el_sep)
    ln_Li_ca = np.log(Li_el_ca)

    # Combine x meshes into single vectors
    x = np.concat([an.x, sep.x, ca.x])
    xm = np.concat([an.xm, sep.xm, ca.xm])
    xp = np.concat([an.xp, sep.xp, ca.xp])

    # Weighted electrolyte properties
    wt_m = 0.5*(xp[:-1] - xm[:-1]) / (x[1:] - x[:-1])
    wt_p = 0.5*(xp[1:] - xm[1:]) / (x[1:] - x[:-1])

    D_el = el.get_D(np.concat([Li_el_an, Li_el_sep, Li_el_ca]), T)
    t0 = el.get_t0(np.concat([Li_el_an, Li_el_sep, Li_el_ca]), T)
    gam = el.get_gamma(np.concat([Li_el_an, Li_el_sep, Li_el_ca]), T)
    kap = el.get_kappa(np.concat([Li_el_an, Li_el_sep, Li_el_ca]), T)

    eps_tau = np.concat([
        an.eps_el**an.p_liq * np.ones(an.Nx),
        sep.eps_el**sep.p_liq * np.ones(sep.Nx),
        ca.eps_el**ca.p_liq * np.ones(ca.Nx),
    ])

    D_eff = wt_m*D_el[:-1]*eps_tau[:-1] + wt_p*D_el[1:]*eps_tau[1:]
    t0_b = np.concat([[t0[0]], wt_m*t0[:-1] + wt_p*t0[1:], [t0[-1]]])
    k_eff = wt_m*kap[:-1]*eps_tau[:-1] + wt_p*kap[1:]*eps_tau[1:]
    gam_b = np.concat([[gam[0]], wt_m*gam[:-1] + wt_p*gam[1:], [gam[-1]]])

    # Ionoic (io) current (i) and molar (N) fluxes in electrolyte
    phi_el = np.concat([phi_el_an, phi_el_sep, phi_el_ca])
    Li_el = np.concat([Li_el_an, Li_el_sep, Li_el_ca])
    ln_Li = np.concat([ln_Li_an, ln_Li_sep, ln_Li_ca])

    ip_io = -k_eff*(phi_el[1:] - phi_el[:-1]) / (x[1:] - x[:-1]) \
          - 2 * k_eff*c.R*T / c.F * (1 + gam_b[1:-1]) * (t0_b[1:-1] - 1) \
              * (ln_Li[1:] - ln_Li[:-1]) / (x[1:] - x[:-1])

    Np_io = D_eff*(Li_el[1:] - Li_el[:-1]) / (x[1:] - x[:-1])

    # Reaction terms ----------------------------------------------------------

    # Anode overpotentials and Li+ productions
    eta = phi_an - phi_el_an - (an.get_Eeq(xs_an[:, -1]) + Hyst_an)
    fluxdir_an = -np.sign(eta)

    i0 = an.get_i0(xs_an[:, -1], Li_el_an, T, fluxdir_an)
    sdot_an = i0 / c.F * (  np.exp( an.alpha_a*c.F*eta / c.R / T)
                          - np.exp(-an.alpha_c*c.F*eta / c.R / T)  )

    # Cathode overpotentials and Li+ productions
    eta = phi_ca - phi_el_ca - (ca.get_Eeq(xs_ca[:, -1]) + Hyst_ca)
    fluxdir_ca = -np.sign(eta)

    i0 = ca.get_i0(xs_ca[:, -1], Li_el_ca, T, fluxdir_ca)
    sdot_ca = i0 / c.F * (  np.exp( ca.alpha_a*c.F*eta / c.R / T)
                          - np.exp(-ca.alpha_c*c.F*eta / c.R / T)  )

    # Boundary condition ------------------------------------------------------
    mode = exp['mode']
    units = exp['units']
    value = exp['value']

    # External current [A/m2]
    i_ext = -sdot_ca[-1]*ca.A_s*c.F*(ca.xp[-1] - ca.xm[-1]) \
          - ca.sigma_s*ca.eps_s**ca.p_sol \
              * (phi_ca[-1] - phi_ca[-2]) / (ca.x[-1] - ca.x[-2])

    if mode == 'current' and units == 'A':
        i_ext = value(t) / bat.area
    elif mode == 'current' and units == 'C':
        i_ext = value(t)*bat.cap / bat.area

    voltage_V = phi_ca[-1]
    current_A = i_ext*bat.area
    power_W = current_A*voltage_V

    # Anode -------------------------------------------------------------------

    # Transference numbers for all anode cell boundaries
    t0_an = t0_b[:an.Nx + 1]

    # Electrolyte fluxes with boundary conditions at x = 0
    im_el = np.concat([[0.], ip_io[:an.Nx - 1]])
    Nm_el = np.concat([[0.], Np_io[:an.Nx - 1]])

    ip_el = ip_io[:an.Nx]
    Np_el = Np_io[:an.Nx]

    # Solid-phase currents w/ BCs at x = 0 and x = an.thick
    s_eff = an.sigma_s*an.eps_s**an.p_sol
    ip_ed = -s_eff*grad_x(an.x, phi_an)

    im_ed = np.concat([[i_ext], ip_ed])
    ip_ed = np.concat([ip_ed, [0.]])

    # Weighted solid particle properties
    wt_m = 0.5*(an.rp[:-1] - an.rm[:-1]) / (an.r[1:] - an.r[:-1])
    wt_p = 0.5*(an.rp[1:] - an.rm[1:]) / (an.r[1:] - an.r[:-1])

    Ds = wt_m*an.get_Ds(xs_an[:, :-1], T, fluxdir_an) \
       + wt_p*an.get_Ds(xs_an[:, 1:], T, fluxdir_an)

    # Solid-phase radial diffusion
    Nk_ed = np.column_stack([
        np.zeros(an.Nx), Ds*grad_r(an.r, Li_an), -sdot_an,
    ])

    fk_ode = div_r(an.rm, an.rp, Nk_ed)

    xr_ptr = an.xr_ptr['Li_ed'].flatten()
    res[xr_ptr] = an.Li_max*svdot[xr_ptr] - fk_ode.flatten()

    # Solid-phase COC (algebraic)
    res[an.x_ptr['phi_ed']] = (ip_ed - im_ed) / (an.xp - an.xm) \
                            + an.A_s*sdot_an*c.F

    # Reference potential BC (algebraic)
    res[an.x_ptr['phi_ed'][0]] = phi_an[0] - 0.

    # Electrolyte COM (differential)
    res[an.x_ptr['Li_el']] = an.eps_el*svdot[an.x_ptr['Li_el']] \
        - (Np_el - Nm_el - (ip_el*t0_an[1:] - im_el*t0_an[:-1]) / c.F) \
        / (an.xp - an.xm) - an.A_s*sdot_an

    # Electrolyte COC (algebraic)
    res[an.x_ptr['phi_el']] = (ip_el - im_el) / (an.xp - an.xm) \
                            - an.A_s*sdot_an*c.F

    # Hysteresis (differential)
    if 'Hysteresis' in an._submodels:
        res[an.x_ptr['hyst']] = svdot[an.x_ptr['hyst']] \
            - np.abs(sdot_an*c.F*an.g_hyst / 3600. / bat.cap) \
            * (sign(sdot_an) - hyst_an)

    # Store some outputs for verification
    if mode == 'post':
        div_i_an = (ip_ed - im_ed + ip_el - im_el) / (an.xp - an.xm)
        sum_ip = ip_el + ip_ed
        i_el_x = im_el

    # Separator ---------------------------------------------------------------

    # Transference numbers for all separator cell boundaries
    t0_sep = t0_b[an.Nx:an.Nx + sep.Nx + 1]

    # Electrolyte fluxes
    im_el = ip_io[an.Nx - 1:an.Nx + sep.Nx - 1]
    Nm_el = Np_io[an.Nx - 1:an.Nx + sep.Nx - 1]

    ip_el = ip_io[an.Nx:an.Nx + sep.Nx]
    Np_el = Np_io[an.Nx:an.Nx + sep.Nx]

    # Electrolyte COC (algebraic)
    res[sep.x_ptr['phi_el']] = (ip_el - im_el) / (sep.xp - sep.xm)

    # Electrolyte COM (differential)
    res[sep.x_ptr['Li_el']] = sep.eps_el*svdot[sep.x_ptr['Li_el']] \
        - (Np_el - Nm_el - (ip_el*t0_sep[1:] - im_el*t0_sep[:-1]) / c.F) \
        / (sep.xp - sep.xm)

    # Store some outputs for verification
    if mode == 'post':
        div_i_sep = (ip_el - im_el) / (sep.xp - sep.xm)
        sum_ip = np.concat([sum_ip, ip_el])
        i_el_x = np.concat([i_el_x, im_el])

    # Cathode -----------------------------------------------------------------

    # Transference numbers for all cathode cell boundaries
    t0_ca = t0_b[an.Nx + sep.Nx:]

    # Electrolyte fluxes with boundary conditions at x = ca.thick
    im_el = ip_io[an.Nx + sep.Nx - 1:an.Nx + sep.Nx + ca.Nx - 1]
    Nm_el = Np_io[an.Nx + sep.Nx - 1:an.Nx + sep.Nx + ca.Nx - 1]

    ip_el = np.concat([ip_io[an.Nx + sep.Nx:], [0.]])
    Np_el = np.concat([Np_io[an.Nx + sep.Nx:], [0.]])

    # Solid-phase currents w/ BCs at x = sep.thick and x = ca.thick
    s_eff = ca.sigma_s*ca.eps_s**ca.p_sol
    ip_ed = -s_eff*grad_x(ca.x, phi_ca)

    im_ed = np.concat([[0.], ip_ed])
    ip_ed = np.concat([ip_ed, [i_ext]])

    # Weighted solid particle properties
    wt_m = 0.5*(ca.rp[:-1] - ca.rm[:-1]) / (ca.r[1:] - ca.r[:-1])
    wt_p = 0.5*(ca.rp[1:] - ca.rm[1:]) / (ca.r[1:] - ca.r[:-1])

    Ds = wt_m*ca.get_Ds(xs_ca[:, :-1], T, fluxdir_ca) \
       + wt_p*ca.get_Ds(xs_ca[:, 1:], T, fluxdir_ca)

    # Solid-phase radial diffusion
    Nk_ed = np.column_stack([
        np.zeros(ca.Nx), Ds*grad_r(ca.r, Li_ca), -sdot_ca,
    ])

    fk_ode = div_r(ca.rm, ca.rp, Nk_ed)

    xr_ptr = ca.xr_ptr['Li_ed'].flatten()
    res[xr_ptr] = ca.Li_max*svdot[xr_ptr] - fk_ode.flatten()

    # Solid-phase COC (algebraic)
    res[ca.x_ptr['phi_ed']] = (ip_ed - im_ed) / (ca.xp - ca.xm) \
                            + ca.A_s*sdot_ca*c.F

    if mode == 'voltage':
        res[ca.x_ptr['phi_ed'][-1]] = voltage_V - value(t)
    elif mode == 'power':
        res[ca.x_ptr['phi_ed'][-1]] = power_W - value(t)

    # Electrolyte COM (differential)
    res[ca.x_ptr['Li_el']] = ca.eps_el*svdot[ca.x_ptr['Li_el']] \
        - (Np_el - Nm_el - (ip_el*t0_ca[1:] - im_el*t0_ca[:-1]) / c.F) \
        / (ca.xp - ca.xm) - ca.A_s*sdot_ca

    # Electrolyte COC (algebraic)
    res[ca.x_ptr['phi_el']] = (ip_el - im_el) / (ca.xp - ca.xm) \
                            - ca.A_s*sdot_ca*c.F

    # Hysteresis (differential)
    if 'Hysteresis' in ca._submodels:
        res[ca.x_ptr['hyst']] = svdot[ca.x_ptr['hyst']] \
            - np.abs(sdot_ca*c.F*ca.g_hyst / 3600. / bat.cap) \
            * (sign(sdot_ca) - hyst_ca)

    # Store some outputs for verification
    if mode == 'post':
        div_i_ca = (ip_ed - im_ed + ip_el - im_el) / (ca.xp - ca.xm)
        sum_ip = np.concat([sum_ip, ip_el + ip_ed])
        i_el_x = np.concat([i_el_x, im_el, [ip_el[-1]]])

    # Events tracking ---------------------------------------------------------
    total_time = sim._t0 + t

    exp['events'] = {
        'time_s': total_time,
        'time_min': total_time / 60.,
        'time_h': total_time / 3600.,
        'current_A': current_A,
        'current_C': current_A / bat.cap,
        'voltage_V': voltage_V,
        'power_W': power_W,
    }

    # Returns -----------------------------------------------------------------
    if mode == 'post':
        return div_i_an, div_i_sep, div_i_ca, sdot_an, sdot_ca, sum_ip, i_el_x
