"""Unit tests for the WindWaves component."""

from nose.tools import assert_is_instance, assert_equal
from oceanwaves.components.windwaves import WindWaves


COMPONENT_NAME = 'WindWaves'


def test_instantiate():
    """Test instantiating WindWaves."""
    x = WindWaves()
    assert_is_instance(x, WindWaves)


def test_finalize():
    """Test BMI finalize method."""
    WindWaves().finalize()


def test_get_component_name():
    """Test BMI get_component_name method."""
    assert_equal(COMPONENT_NAME, WindWaves().get_component_name())
