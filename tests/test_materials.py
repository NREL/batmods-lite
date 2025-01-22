import os

import pytest
import numpy as np
import bmlite as bm


@pytest.fixture(scope='module')
def args():
    alpha_a = 0.5
    alpha_c = 0.5
    Li_max = 50.
    return alpha_a, alpha_c, Li_max


def test_gen2_electrolyte():
    el = bm.materials.Gen2Electrolyte()

    directory = os.path.dirname(__file__)
    data = np.load(directory + '/materials_data/gen2electrolyte.npz')

    C_Li, T = data['C_Li'], data['T']

    D = el.get_D(C_Li, T)
    assert np.allclose(D / max(D), data['D'] / max(D))

    t0 = el.get_t0(C_Li, T)
    assert np.allclose(t0 / max(t0), data['t0'] / max(t0))

    kappa = el.get_kappa(C_Li, T)
    assert np.allclose(kappa / max(kappa), data['kappa'] / max(kappa))

    gamma = el.get_gamma(C_Li, T)
    assert np.allclose(gamma / max(gamma), data['gamma'] / max(gamma))

    data.close()


def test_graphite_fast(args):
    gr = bm.materials.GraphiteFast(args[0], args[1], args[2])

    directory = os.path.dirname(__file__)
    data = np.load(directory + '/materials_data/graphitefast.npz')

    x, C_Li, T = data['x'], data['C_Li'], data['T']

    Ds = gr.get_Ds(x, T)
    assert np.allclose(Ds / max(Ds), data['Ds'] / max(Ds))

    i0 = gr.get_i0(x, C_Li, T)
    assert np.allclose(i0 / max(i0), data['i0'] / max(i0))

    Eeq = gr.get_Eeq(x)
    assert np.allclose(Eeq / max(Eeq), data['Eeq'] / max(Eeq))

    data.close()


def test_graphite_slow(args):
    gr = bm.materials.GraphiteSlow(args[0], args[1], args[2])

    directory = os.path.dirname(__file__)
    data = np.load(directory + '/materials_data/graphiteslow.npz')

    x, C_Li, T = data['x'], data['C_Li'], data['T']

    Ds = gr.get_Ds(x, T)
    assert np.allclose(Ds / max(Ds), data['Ds'] / max(Ds))

    i0 = gr.get_i0(x, C_Li, T)
    assert np.allclose(i0 / max(i0), data['i0'] / max(i0))

    Eeq = gr.get_Eeq(x)
    assert np.allclose(Eeq / max(Eeq), data['Eeq'] / max(Eeq))

    data.close()

    with pytest.raises(ValueError):
        gr.get_Eeq(gr.x_min - 0.01)

    with pytest.raises(ValueError):
        gr.get_Eeq(gr.x_max * np.ones(5) + 0.01)


def test_nmc_532_fast(args):
    nmc = bm.materials.NMC532Fast(args[0], args[1], args[2])

    directory = os.path.dirname(__file__)
    data = np.load(directory + '/materials_data/nmc532fast.npz')

    x, C_Li, T = data['x'], data['C_Li'], data['T']

    Ds = nmc.get_Ds(x, T)
    assert np.allclose(Ds / max(Ds), data['Ds'] / max(Ds))

    i0 = nmc.get_i0(x, C_Li, T)
    assert np.allclose(i0 / max(i0), data['i0'] / max(i0))

    Eeq = nmc.get_Eeq(x)
    assert np.allclose(Eeq / max(Eeq), data['Eeq'] / max(Eeq))

    data.close()


def test_nmc_532_slow(args):
    nmc = bm.materials.NMC532Slow(args[0], args[1], args[2])

    directory = os.path.dirname(__file__)
    data = np.load(directory + '/materials_data/nmc532slow.npz')

    x, C_Li, T = data['x'], data['C_Li'], data['T']

    Ds = nmc.get_Ds(x, T)
    assert np.allclose(Ds / max(Ds), data['Ds'] / max(Ds))

    i0 = nmc.get_i0(x, C_Li, T)
    assert np.allclose(i0 / max(i0), data['i0'] / max(i0))

    Eeq = nmc.get_Eeq(x)
    assert np.allclose(Eeq / max(Eeq), data['Eeq'] / max(Eeq))

    data.close()

    with pytest.raises(ValueError):
        nmc.get_Eeq(nmc.x_min - 0.01)

    with pytest.raises(ValueError):
        nmc.get_Eeq(nmc.x_max * np.ones(5) + 0.01)
