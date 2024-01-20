import pytest

import bmlite as bm


@pytest.fixture(scope='session')
def sim():
    sim = bm.SPM.Simulation()
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
    bm.SPM.templates()
    bm.SPM.templates(0)
    bm.SPM.templates('default_SPM')
    bm.SPM.templates('default_SPM.yaml')
    bm.SPM.templates(exp=0)
    bm.SPM.templates(exp='constant_current')
    bm.SPM.templates(exp='constant_current.yaml')
    assert True
