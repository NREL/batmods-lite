import pytest

import numpy as np
import bmlite as bm


class XDomain():

    # Mimics a P2D separator
    def __init__(self):
        self.Nx = 10

        self.ptr = {}
        self.ptr['phie'] = 0
        self.ptr['ce'] = 1

        self.ptr['x_off'] = 2


class RDomain():

    # Mimics a SPM electrode
    def __init__(self):
        self.Nr = 10

        self.ptr = {}
        self.ptr['xs'] = 0
        self.ptr['phis'] = self.Nr

        self.ptr['r_off'] = 1


class XRDomain():

    # Mimics a P2D electrode
    def __init__(self):
        self.Nr = 10
        self.Nx = 10

        self.ptr = {}
        self.ptr['xs'] = 0
        self.ptr['phis'] = self.Nr
        self.ptr['ce'] = self.ptr['phis'] + 1
        self.ptr['phie'] = self.ptr['ce'] + 1

        self.ptr['r_off'] = 1
        self.ptr['x_off'] = self.Nr + 3


@pytest.fixture(scope='function')
def x_domain():
    return XDomain()


@pytest.fixture(scope='function')
def r_domain():
    return RDomain()


@pytest.fixture(scope='function')
def xr_domain():
    return XRDomain()


def test_x_ptr(x_domain):

    bm.mesh.x_ptr(x_domain, ['ce', 'phie'])

    ce_ptr = np.array([x_domain.ptr['ce'] + i * x_domain.ptr['x_off']
                       for i in range(x_domain.Nx)], dtype=int)

    phie_ptr = np.array([x_domain.ptr['phie'] + i * x_domain.ptr['x_off']
                         for i in range(x_domain.Nx)], dtype=int)

    assert np.all(x_domain.x_ptr['ce'] == ce_ptr)
    assert np.all(x_domain.x_ptr['phie'] == phie_ptr)


def test_r_ptr(r_domain):

    bm.mesh.r_ptr(r_domain, ['xs'])

    cs_ptr = np.array([r_domain.ptr['xs'] + j * r_domain.ptr['r_off']
                       for j in range(r_domain.Nr)], dtype=int)

    assert np.all(r_domain.r_ptr['xs'] == cs_ptr)


def test_xr_ptr(xr_domain):

    bm.mesh.xr_ptr(xr_domain, ['xs'])

    cs_ptr = np.zeros([xr_domain.Nx, xr_domain.Nr], dtype=int)
    for i in range(xr_domain.Nx):
        for j in range(xr_domain.Nr):
            cs_ptr[i, j] = xr_domain.ptr['xs'] + i * xr_domain.ptr['x_off'] \
                + j * xr_domain.ptr['r_off']

    assert np.all(xr_domain.xr_ptr['xs'] == cs_ptr)


def test_uniform_mesh():

    xm, xp, x = bm.mesh.uniform_mesh(10., 50, 5.)

    dx = 10. / 50

    assert np.allclose(xm, np.linspace(5., 5. + 10. - dx, 50))
    assert np.allclose(xp, np.linspace(5. + dx, 5. + 10., 50))
    assert np.allclose(x, 0.5 * (np.linspace(5., 5. + 10. - dx, 50)
                                 + np.linspace(5. + dx, 5. + 10., 50)))


def test_param_weights():

    xm = np.linspace(0, 5, 6)**2
    xp = np.linspace(1, 6, 6)**2

    wt_m1, wt_p1 = bm.mesh.param_weights(xm, xp)

    x = 0.5 * (xm + xp)

    wt_m2 = 0.5 * (xp[: -1] - xm[: -1]) / (x[1:] - x[: -1])
    wt_p2 = 0.5 * (xp[1:] - xm[1:]) / (x[1:] - x[: -1])

    assert np.allclose(wt_m1, wt_m2)
    assert np.allclose(wt_p1, wt_p2)
