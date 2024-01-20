"""
DAE Module
----------
This module includes the system of differential algebraic equations (DAE) for
the SPM model. In addition, the ``bandwidth`` function is included in this
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
    inputs : SPM Simulation object
        An instance of the SPM model simulation. See
        :class:`bmlite.SPM.Simulation`.

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
    exp['C_rate'] = 0.
    exp['t_min'] = 0.
    exp['t_max'] = 3600.
    exp['Nt'] = 3601

    exp['BC'] = 'current'
    exp['i_ext'] = exp['C_rate'] * sim.bat.cap / sim.bat.area

    # Turn on band flag and set BC
    sim._flags['band'] = True
    sim._flags['BC'] = 'current'

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
        dres = res_0 - residuals(0., sv, dsvp, res, (sim, exp))

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

    # Turn off band flag and reset BC
    sim._flags['band'] = False
    sim._flags['BC'] = None

    # Make Jacobian pattern of zeros and ones
    j_pat = np.zeros_like(jac)
    j_pat[jac != 0] = 1

    return lband, uband, j_pat


def residuals(t: float, sv: _ndarray, svdot: _ndarray, res: _ndarray,
              inputs: tuple[object, dict]) -> None | tuple[_ndarray]:
    """
    The DAE residuals ``res = M*y' - f(t, y)`` for the SPM model.

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

    inputs : (sim : SPM Simulation object, exp : experiment dict)
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

        ========== =======================================================
        Variable   Description [units] (type)
        ========== =======================================================
        res        residuals ``res = M*y' - f(t, y)`` [units] (*1D array*)
        sdot_an    anode Li+ production rate [kmol/m^3/s] (*float*)
        sdot_ca    cathode Li+ production rate [kmol/m^3/s] (*float*)
        ========== =======================================================
    """

    import numpy as np

    from .. import Constants

    c = Constants()

    # Break inputs into separate objects
    sim, exp = inputs

    bat, el, an, ca = sim.bat, sim.el, sim.an, sim.ca

    # Simulation temperature
    T = bat.temp

    # Organize values from sv
    phi_an = sv[an.ptr['phi_ed']]
    phi_el = sv[el.ptr['phi_el']]
    phi_ca = sv[ca.ptr['phi_ed']]

    Li_an = sv[an.r_ptr('Li_ed')] * an.Li_max
    Li_ca = sv[ca.r_ptr('Li_ed')] * ca.Li_max

    # Anode -------------------------------------------------------------------

    # Reaction current
    eta = phi_an - phi_el - an.get_Eeq(Li_an[-1] / an.Li_max, T)

    i0 = an.get_i0(Li_an[-1] / an.Li_max, el.Li_0, T)
    sdot_an = i0 / c.F * (  np.exp( an.alpha_a * c.F * eta / c.R / T)
                          - np.exp(-an.alpha_c * c.F * eta / c.R / T)  )

    # Weighted solid particle properties
    wt_m = 0.5 * (an.rp[:-1] - an.rm[:-1]) / (an.r[1:] - an.r[:-1])
    wt_p = 0.5 * (an.rp[1:] - an.rm[1:]) / (an.r[1:] - an.r[:-1])

    Ds_an = wt_m * an.get_Ds(Li_an[:-1] / an.Li_max, T) \
          + wt_p * an.get_Ds(Li_an[1:] / an.Li_max, T)

    # Solid-phase COM (differential)
    Nm_ed = np.zeros(an.Nr)
    Np_ed = np.zeros(an.Nr)

    Np_ed[:-1] = Ds_an * (Li_an[1:] - Li_an[:-1]) / (an.r[1:] - an.r[:-1])
    Nm_ed[1:] = Np_ed[:-1]

    Np_ed[-1] = -sdot_an

    res[an.r_ptr('Li_ed')] = an.Li_max * svdot[an.r_ptr('Li_ed')] \
                           - ( an.rp**2 * Np_ed - an.rm**2 * Nm_ed ) \
                             / an.r**2 / (an.rp - an.rm)

    # Solid-phase COC (algebraic)
    res[an.ptr['phi_ed']] = phi_an - 0.

    # Cathode -----------------------------------------------------------------

    # Reaction current
    eta = phi_ca - phi_el - ca.get_Eeq(Li_ca[-1] / ca.Li_max, T)

    i0 = ca.get_i0(Li_ca[-1] / ca.Li_max, el.Li_0, T)
    sdot_ca = i0 / c.F * (  np.exp( ca.alpha_a * c.F * eta / c.R / T)
                          - np.exp(-ca.alpha_c * c.F * eta / c.R / T)  )

    # Weighted solid particle properties
    wt_m = 0.5 * (ca.rp[:-1] - ca.rm[:-1]) / (ca.r[1:] - ca.r[:-1])
    wt_p = 0.5 * (ca.rp[1:] - ca.rm[1:]) / (ca.r[1:] - ca.r[:-1])

    Ds_ca = wt_m * ca.get_Ds(Li_ca[:-1] / ca.Li_max, T) \
          + wt_p * ca.get_Ds(Li_ca[1:] / ca.Li_max, T)

    # Solid-phase COM (differential)
    Nm_ed = np.zeros(ca.Nr)
    Np_ed = np.zeros(ca.Nr)

    Np_ed[:-1] = Ds_ca * (Li_ca[1:] - Li_ca[:-1]) / (ca.r[1:] - ca.r[:-1])
    Nm_ed[1:] = Np_ed[:-1]

    Np_ed[-1] = -sdot_ca

    res[ca.r_ptr('Li_ed')] = ca.Li_max * svdot[ca.r_ptr('Li_ed')] \
                           - ( ca.rp**2 * Np_ed - ca.rm**2 * Nm_ed ) \
                             / ca.r**2 / (ca.rp - ca.rm)

    # External current [A/m^2]
    exp['i_ext'] = -sdot_an * an.A_s * an.thick * c.F

    # Solid-phase COC (algebraic)
    if sim._flags['BC'] == 'current':
        res[ca.ptr['phi_ed']] = sdot_ca * ca.A_s * ca.thick * c.F \
                              - bat.cap * exp['C_rate'] / bat.area

    elif sim._flags['BC'] == 'voltage':
        res[ca.ptr['phi_ed']] = phi_ca - exp['V_ext']

    elif sim._flags['BC'] == 'power':
        res[ca.ptr['phi_ed']] = exp['i_ext'] * phi_ca - exp['P_ext']

    # Electrolyte -------------------------------------------------------------

    # Electrolyte potential (algebraic)
    if sim._flags['BC'] == 'current':
        res[el.ptr['phi_el']] = sdot_an * an.A_s * an.thick * c.F \
                              + bat.cap * exp['C_rate'] / bat.area

    elif sim._flags['BC'] == 'voltage':
        res[el.ptr['phi_el']] = sdot_an * an.A_s * an.thick \
                              + sdot_ca * ca.A_s * ca.thick

    elif sim._flags['BC'] == 'power':
        res[el.ptr['phi_el']] = sdot_an * an.A_s * an.thick \
                              + sdot_ca * ca.A_s * ca.thick

    # Returns -----------------------------------------------------------------

    if sim._flags['band']:
        return res

    elif sim._flags['post']:
        return res, sdot_an, sdot_ca
