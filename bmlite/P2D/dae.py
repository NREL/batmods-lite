"""
DAE Module
----------
This module includes the system of differential algebraic equations (DAE) for
the P2D model. In addition, the ``bandwidth`` function is included in this
module, which helps determine the lower and upper bandwidths of ``residuals``
so the ``'band'`` linear solver option can be used in the ``IDASolver`` class.
"""

from numpy import ndarray as _ndarray


def bandwidth(sim: object) -> tuple[int | _ndarray]:
    """
    Determine the DAE system's bandwidth and Jacobian pattern.

    Numerically determines the bandwidth and Jacobian pattern of the residual
    function by perturbating each ``y`` and ``yp`` term and determining which
    ``dres/dy`` and ``dres/dy'`` terms are non-zero. The bandwidth is required
    to use the "band" option for IDA's linear solver, which speeds up each
    integration step compared to the "dense" linear solver.

    Parameters
    ----------
    inputs : P2D Simulation object
        An instance of the P2D model simulation. See
        :class:`bmlite.P2D.Simulation`.

    Returns
    -------
    lband : int
        Lower bandwidth from the residual function's Jacobian pattern.

    uband : int
        Upper bandwidth from the residual function's Jacobian pattern.

    j_pat : 2D array
        Residual function Jacobian pattern, as an array of ones and zeros.
    """

    import numpy as np

    # Jacobian size
    N = sim.sv_0.size

    # Fake OCV experiment
    exp = {}
    exp["C_rate"] = 0.0
    exp["t_min"] = 0.0
    exp["t_max"] = 3600.0
    exp["Nt"] = 3601

    exp["BC"] = "current"
    exp["i_ext"] = exp["C_rate"] * sim.bat.cap / sim.bat.area

    # Turn on band flag
    sim._flags["band"] = True

    # Perturbed variables
    jac = np.zeros([N, N])
    sv = sim.sv_0.copy()
    svdot = sim.svdot_0.copy()
    res = np.zeros_like(sv)
    res_0 = residuals(0, sv, svdot, res, (sim, exp))

    for j in range(N):
        dsv = np.copy(sv)
        res = np.copy(res)
        svdot = np.copy(svdot)

        dsv[j] = 1.01 * (sv[j] + 0.01)
        dres = res_0 - residuals(0, dsv, svdot, res, (sim, exp))

        jac[:, j] = dres

    for j in range(N):
        sv = np.copy(sv)
        res = np.copy(res)
        dsvp = np.copy(svdot)

        dsvp[j] = 1.01 * (svdot[j] + 0.01)
        dres = res_0 - residuals(0.0, sv, dsvp, res, (sim, exp))

        jac[:, j] += dres

    # Find lband and uband
    lband = 0
    uband = 0

    for i in range(jac.shape[0]):
        l_inds = np.where(abs(jac[i, :i]) > 0)[0]
        if len(l_inds) >= 1 and i - l_inds[0] > lband:
            lband = i - l_inds[0]

        u_inds = i + np.where(abs(jac[i, i:]) > 0)[0]
        if len(u_inds) >= 1 and u_inds[-1] - i > uband:
            uband = u_inds[-1] - i

    # Turn off band flag
    sim._flags["band"] = False

    # Make Jacobian pattern of zeros and ones
    j_pat = np.zeros_like(jac)
    j_pat[jac != 0] = 1

    return lband, uband, j_pat


def residuals(
    t: float,
    sv: _ndarray,
    svdot: _ndarray,
    res: _ndarray,
    inputs: tuple[object, dict],
) -> None | tuple[_ndarray]:
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
    None
        If no ``sim._flags`` are ``True``.

    res : 1D array
        Array of residuals if ``sim._flags['band'] = True``.

    outputs : tuple[1D array]
        If ``sim._flags['post'] = True`` then ``outputs`` is returned, which
        includes post-processed values. These can help verify the governing
        equations and boundary conditions are satisfied. They can also be
        useful for interpreting causes of good/bad battery performance. The
        order and description of the arrays is given below:

        ======== =============================================================
        Variable Description [units] (type)
        ======== =============================================================
        res      residuals ``res = M*y' - f(t, y)`` [units] (*1D array*)
        sdot_an  Li+ production rate at each x_a [kmol/m^3/s] (*1D array*)
        sdot_ca  Li+ production rate at each x_c [kmol/m^3/s] (*1D array*)
        sum_ip   ``i_ed + i_el`` at each "plus" interface [A/m^2] (*1D array*)
        i_el_x   ``i_el`` at each interface [A/m^2] (*1D array*)
        ======== =============================================================
    """

    import numpy as np

    from .. import Constants

    c = Constants()

    # Break inputs into separate objects
    sim, exp = inputs

    bat, el, an, sep, ca = sim.bat, sim.el, sim.an, sim.sep, sim.ca

    # Simulation temperature
    T = bat.temp

    # Organize values from sv
    Li_an = np.zeros([an.Nx, an.Nr])
    phi_an = sv[an.x_ptr("phi_ed")]
    Li_el_an = sv[an.x_ptr("Li_el")]
    phi_el_an = sv[an.x_ptr("phi_el")]

    Li_el_sep = sv[sep.x_ptr("Li_el")]
    phi_el_sep = sv[sep.x_ptr("phi_el")]

    Li_ca = np.zeros([ca.Nx, ca.Nr])
    phi_ca = sv[ca.x_ptr("phi_ed")]
    Li_el_ca = sv[ca.x_ptr("Li_el")]
    phi_el_ca = sv[ca.x_ptr("phi_el")]

    # sv values for anode particles (both x and r)
    for j in range(an.Nr):
        r_off = j * an.ptr["r_off"]
        Li_an[:, j] = sv[an.x_ptr("Li_ed", r_off)] * an.Li_max

    # sv values for cathode particles (both x and r)
    for j in range(ca.Nr):
        r_off = j * ca.ptr["r_off"]
        Li_ca[:, j] = sv[ca.x_ptr("Li_ed", r_off)] * ca.Li_max

    # Pre-calculate ln(Li_el)
    ln_Li_an = np.log(Li_el_an)
    ln_Li_sep = np.log(Li_el_sep)
    ln_Li_ca = np.log(Li_el_ca)

    # Combine x meshes into single vectors
    x_all = np.hstack([an.x, sep.x, ca.x])
    xm_all = np.hstack([an.xm, sep.xm, ca.xm])
    xp_all = np.hstack([an.xp, sep.xp, ca.xp])

    # Weighted electrolyte properties
    wt_m = 0.5 * (xp_all[:-1] - xm_all[:-1]) / (x_all[1:] - x_all[:-1])
    wt_p = 0.5 * (xp_all[1:] - xm_all[1:]) / (x_all[1:] - x_all[:-1])

    D_el = el.get_D(np.hstack([Li_el_an, Li_el_sep, Li_el_ca]), T)
    t0 = el.get_t0(np.hstack([Li_el_an, Li_el_sep, Li_el_ca]), T)
    gam = el.get_gamma(np.hstack([Li_el_an, Li_el_sep, Li_el_ca]), T)
    kap = el.get_kappa(np.hstack([Li_el_an, Li_el_sep, Li_el_ca]), T)

    eps_tau = np.hstack(
        [
            an.eps_el**an.p_liq * np.ones(an.Nx),
            sep.eps_el**sep.p_liq * np.ones(sep.Nx),
            ca.eps_el**ca.p_liq * np.ones(ca.Nx),
        ]
    )

    D_eff = wt_m * D_el[:-1] * eps_tau[:-1] + wt_p * D_el[1:] * eps_tau[1:]
    t0_b = np.hstack([t0[0], wt_m * t0[:-1] + wt_p * t0[1:], t0[-1]])
    k_eff = wt_m * kap[:-1] * eps_tau[:-1] + wt_p * kap[1:] * eps_tau[1:]
    gam_b = np.hstack([gam[0], wt_m * gam[:-1] + wt_p * gam[1:], gam[-1]])

    # Ionoic (io) current (i) and molar (N) fluxes in electrolyte
    phi_el = np.hstack([phi_el_an, phi_el_sep, phi_el_ca])
    Li_el = np.hstack([Li_el_an, Li_el_sep, Li_el_ca])
    ln_Li = np.hstack([ln_Li_an, ln_Li_sep, ln_Li_ca])

    ip_io = -k_eff * (phi_el[1:] - phi_el[:-1]) / (
        x_all[1:] - x_all[:-1]
    ) - 2 * k_eff * c.R * T / c.F * (1 + gam_b[1:-1]) * (t0_b[1:-1] - 1) * (
        ln_Li[1:] - ln_Li[:-1]
    ) / (
        x_all[1:] - x_all[:-1]
    )

    Np_io = D_eff * (Li_el[1:] - Li_el[:-1]) / (x_all[1:] - x_all[:-1])

    # Anode--------------------------------------------------------------------

    # Transference numbers for all anode cell boundaries
    t0_an = t0_b[: an.Nx + 1]

    # Electrolyte fluxes with boundary conditions at x = 0
    im_el = np.hstack([0.0, ip_io[: an.Nx - 1]])
    Nm_el = np.hstack([0.0, Np_io[: an.Nx - 1]])

    ip_el = ip_io[: an.Nx]
    Np_el = Np_io[: an.Nx]

    # Solid-phase currents w/ BCs at x = 0 and x = an.thick
    s_eff = an.sigma_s * an.eps_s**an.p_sol
    ip_ed = -s_eff * (phi_an[1:] - phi_an[:-1]) / (an.x[1:] - an.x[:-1])

    im_ed = np.hstack([-exp["i_ext"], ip_ed])
    ip_ed = np.hstack([ip_ed, 0.0])

    # Overpotentials and Li+ productions
    eta = phi_an - phi_el_an - an.get_Eeq(Li_an[:, -1] / an.Li_max, T)

    i0 = an.get_i0(Li_an[:, -1] / an.Li_max, Li_el_an, T)
    sdot_an = (
        i0
        * (
            np.exp(an.alpha_a * c.F * eta / c.R / T)
            - np.exp(-an.alpha_c * c.F * eta / c.R / T)
        )
        / c.F
    )

    # Weighted solid particle properties
    wt_m = 0.5 * (an.rp[:-1] - an.rm[:-1]) / (an.r[1:] - an.r[:-1])
    wt_p = 0.5 * (an.rp[1:] - an.rm[1:]) / (an.r[1:] - an.r[:-1])

    Ds = wt_m * an.get_Ds(Li_an[:, :-1] / an.Li_max, T) + wt_p * an.get_Ds(
        Li_an[:, 1:] / an.Li_max, T
    )

    # Solid-phase radial diffusion
    for j in range(an.Nr):
        if j == 0:
            Nm_ed = 0.0
            Np_ed = (
                Ds[:, j]
                * (Li_an[:, j + 1] - Li_an[:, j])
                / (an.r[j + 1] - an.r[j])
            )

        elif j < an.Nr - 1:
            Nm_ed = Np_ed
            Np_ed = (
                Ds[:, j]
                * (Li_an[:, j + 1] - Li_an[:, j])
                / (an.r[j + 1] - an.r[j])
            )

        elif j == an.Nr - 1:
            Nm_ed = Np_ed
            Np_ed = -sdot_an

        # Solid-phase COM (differential)
        ptr = an.x_ptr("Li_ed", j * an.ptr["r_off"])
        res[ptr] = an.Li_max * svdot[ptr] - 1 / an.r[j] ** 2 * 1 / (
            an.rp[j] - an.rm[j]
        ) * (an.rp[j] ** 2 * Np_ed - an.rm[j] ** 2 * Nm_ed)

    # Solid-phase COC (algebraic)
    res[an.x_ptr("phi_ed")] = (ip_ed - im_ed) / (
        an.xp - an.xm
    ) + an.A_s * sdot_an * c.F

    # Reference potential BC (algebraic)
    res[an.x_ptr("phi_ed")[0]] = phi_an[0]

    # Electrolyte COM (differential)
    res[an.x_ptr("Li_el")] = (
        an.eps_el * svdot[an.x_ptr("Li_el")]
        - (Np_el - Nm_el - (ip_el * t0_an[1:] - im_el * t0_an[:-1]) / c.F)
        / (an.xp - an.xm)
        - an.A_s * sdot_an
    )

    # Electrolyte COC (algebraic)
    res[an.x_ptr("phi_el")] = (ip_el - im_el) / (
        an.xp - an.xm
    ) - an.A_s * sdot_an * c.F

    # Store some outputs for verification
    sum_ip = ip_el + ip_ed
    i_el_x = im_el

    # Separator----------------------------------------------------------------

    # Transference numbers for all separator cell boundaries
    t0_sep = t0_b[an.Nx : an.Nx + sep.Nx + 1]

    # Electrolyte fluxes
    im_el = ip_io[an.Nx - 1 : an.Nx + sep.Nx - 1]
    Nm_el = Np_io[an.Nx - 1 : an.Nx + sep.Nx - 1]

    ip_el = ip_io[an.Nx : an.Nx + sep.Nx]
    Np_el = Np_io[an.Nx : an.Nx + sep.Nx]

    # Electrolyte COC (algebraic)
    res[sep.x_ptr("phi_el")] = (ip_el - im_el) / (sep.xp - sep.xm)

    # Electrolyte COM (differential)
    res[sep.x_ptr("Li_el")] = sep.eps_el * svdot[sep.x_ptr("Li_el")] - (
        Np_el - Nm_el - (ip_el * t0_sep[1:] - im_el * t0_sep[:-1]) / c.F
    ) / (sep.xp - sep.xm)

    # Store some outputs for verification
    sum_ip = np.hstack([sum_ip, ip_el])
    i_el_x = np.hstack([i_el_x, im_el])

    # Cathode------------------------------------------------------------------

    # Transference numbers for all cathode cell boundaries
    t0_ca = t0_b[an.Nx + sep.Nx :]

    # Electrolyte fluxes with boundary conditions at x = ca.thick
    im_el = ip_io[an.Nx + sep.Nx - 1 : an.Nx + sep.Nx + ca.Nx - 1]
    Nm_el = Np_io[an.Nx + sep.Nx - 1 : an.Nx + sep.Nx + ca.Nx - 1]

    ip_el = np.hstack([ip_io[an.Nx + sep.Nx :], 0.0])
    Np_el = np.hstack([Np_io[an.Nx + sep.Nx :], 0.0])

    # Solid-phase currents w/ BCs at x = sep.thick and x = ca.thick
    s_eff = ca.sigma_s * ca.eps_s**ca.p_sol
    ip_ed = -s_eff * (phi_ca[1:] - phi_ca[:-1]) / (ca.x[1:] - ca.x[:-1])

    im_ed = np.hstack([0.0, ip_ed])
    ip_ed = np.hstack([ip_ed, -exp["i_ext"]])

    # Overpotentials and Li+ productions
    eta = phi_ca - phi_el_ca - ca.get_Eeq(Li_ca[:, -1] / ca.Li_max, T)

    i0 = ca.get_i0(Li_ca[:, -1] / ca.Li_max, Li_el_ca, T)
    sdot_ca = (
        i0
        * (
            np.exp(ca.alpha_a * c.F * eta / c.R / T)
            - np.exp(-ca.alpha_c * c.F * eta / c.R / T)
        )
        / c.F
    )

    # Weighted solid particle properties
    wt_m = 0.5 * (ca.rp[:-1] - ca.rm[:-1]) / (ca.r[1:] - ca.r[:-1])
    wt_p = 0.5 * (ca.rp[1:] - ca.rm[1:]) / (ca.r[1:] - ca.r[:-1])

    Ds = wt_m * ca.get_Ds(Li_ca[:, :-1] / ca.Li_max, T) + wt_p * ca.get_Ds(
        Li_ca[:, 1:] / ca.Li_max, T
    )

    # Solid-phase radial diffusion
    for j in range(ca.Nr):
        if j == 0:
            Nm_ed = 0.0
            Np_ed = (
                Ds[:, j]
                * (Li_ca[:, j + 1] - Li_ca[:, j])
                / (ca.r[j + 1] - ca.r[j])
            )

        elif j < ca.Nr - 1:
            Nm_ed = Np_ed
            Np_ed = (
                Ds[:, j]
                * (Li_ca[:, j + 1] - Li_ca[:, j])
                / (ca.r[j + 1] - ca.r[j])
            )

        elif j == ca.Nr - 1:
            Nm_ed = Np_ed
            Np_ed = -sdot_ca

        # Solid-phase COM (differential)
        ptr = ca.x_ptr("Li_ed", j * ca.ptr["r_off"])
        res[ptr] = ca.Li_max * svdot[ptr] - 1 / ca.r[j] ** 2 * 1 / (
            ca.rp[j] - ca.rm[j]
        ) * (ca.rp[j] ** 2 * Np_ed - ca.rm[j] ** 2 * Nm_ed)

    # Solid-phase COC (algebraic)
    res[ca.x_ptr("phi_ed")] = (ip_ed - im_ed) / (
        ca.xp - ca.xm
    ) + ca.A_s * sdot_ca * c.F

    # Electrolyte COM (differential)
    res[ca.x_ptr("Li_el")] = (
        ca.eps_el * svdot[ca.x_ptr("Li_el")]
        - (Np_el - Nm_el - (ip_el * t0_ca[1:] - im_el * t0_ca[:-1]) / c.F)
        / (ca.xp - ca.xm)
        - ca.A_s * sdot_ca
    )

    # Electrolyte COC (algebraic)
    res[ca.x_ptr("phi_el")] = (ip_el - im_el) / (
        ca.xp - ca.xm
    ) - ca.A_s * sdot_ca * c.F

    # Store some outputs for verification
    sum_ip = np.hstack([sum_ip, ip_el + ip_ed])
    i_el_x = np.hstack([i_el_x, im_el, ip_el[-1]])

    if sim._flags["band"]:
        return res

    elif sim._flags["post"]:
        return res, sdot_an, sdot_ca, sum_ip, i_el_x
