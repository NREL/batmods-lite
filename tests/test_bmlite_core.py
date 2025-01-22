import pytest
import bmlite as bm


@pytest.fixture(scope='module')
def c():
    c = bm.Constants()
    return c


def test_constants_read_only(c):
    with pytest.raises(AttributeError):
        c.F = 1

    with pytest.raises(AttributeError):
        c.R = 1

    with pytest.raises(AttributeError):
        c.a = 1


def test_constants_empty_slots(c):
    assert len(c.__slots__) == 0


def test_constants_faraday(c):
    assert c.F == 96485.3321e3


def test_constants_idealgas(c):
    assert c.R == 8.3145e3


def test_templates():
    with pytest.raises(AttributeError):
        bm.templates('fake')

    bm.templates('spm')
    bm.templates('spm', 0)
    bm.templates('spm', 'graphite_nmc532')
    bm.templates('spm', 'graphite_nmc532.yaml')
    
    bm.templates('p2d')
    bm.templates('p2d', 0)
    bm.templates('p2d', 'graphite_nmc532')
    bm.templates('p2d', 'graphite_nmc532.yaml')

    assert True
