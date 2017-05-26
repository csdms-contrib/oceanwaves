"""Unit tests for theBottomWaveVelocity component. """

from nose.tools import assert_is_instance, assert_equal
from oceanwaves.components.bottomwavevelocity import BottomWaveVelocity


_name = 'BottomWaveVelocity'


def test_instantiate():
    x = BottomWaveVelocity()
    assert_is_instance(x, BottomWaveVelocity)


def test_finalize():
    x = BottomWaveVelocity()
    x.finalize()


def test_get_component_name():
    x = BottomWaveVelocity()
    assert_equal(_name, x.get_component_name())
