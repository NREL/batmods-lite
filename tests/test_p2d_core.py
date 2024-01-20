import pytest

import bmlite as bm


@pytest.fixture(scope='session')
def sim():
    sim = bm.P2D.Simulation()
    return sim


def test_simulation(sim):
    assert sim


def test_j_pattern(sim):
    sim.j_pattern()
    assert True


def test_copy(sim):
    sim2 = sim.copy()
    assert all(sim2.sv_0 == sim.sv_0)


def test_templates():
    bm.P2D.templates()
    bm.P2D.templates(0)
    bm.P2D.templates('default_P2D')
    bm.P2D.templates('default_P2D.yaml')
    bm.P2D.templates(exp=0)
    bm.P2D.templates(exp='constant_current')
    bm.P2D.templates(exp='constant_current.yaml')
    assert True
