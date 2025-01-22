import pytest
import bmlite as bm
import matplotlib.pyplot as plt


@pytest.fixture(scope='module')
def sim():
    sim = bm.SPM.Simulation()
    return sim


def test_simulation(sim):
    assert sim


def test_fake_yaml():
    with pytest.raises(FileNotFoundError):
        _ = bm.SPM.Simulation('fake.yaml')


def test_j_pattern(sim):
    with plt.ioff():
        lband, uband = sim.j_pattern(return_bands=True)

    assert sim._lband == lband
    assert sim._uband == uband

    plt.close('all')


def test_copy(sim):
    sim2 = sim.copy()
    assert all(sim2._sv0 == sim._sv0)
