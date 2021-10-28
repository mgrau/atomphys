import json
import os.path

import pint
import pytest

from atomphys import util
from atomphys.util import default_units, disk_cache, fsolve


def test_default_units():
    class Test:
        _ureg = None
        _quantity = None

        def __init__(self, ureg):
            self._ureg = ureg

        @property
        def test_quantity(self):
            return self._quantity

        @test_quantity.setter
        @default_units("meter")
        def test_quantity(self, new_quantity):
            self._quantity = new_quantity

    ureg = pint.UnitRegistry()
    test_object = Test(ureg)

    test_object.test_quantity = 0
    assert isinstance(test_object.test_quantity, pint.Quantity)
    assert test_object.test_quantity.check("[length]")

    test_object.test_quantity = 1
    assert isinstance(test_object.test_quantity, pint.Quantity)
    assert test_object.test_quantity.check("[length]")

    test_object.test_quantity = "1 m"
    assert isinstance(test_object.test_quantity, pint.Quantity)
    assert test_object.test_quantity.check("[length]")
    assert test_object.test_quantity.m == 1

    test_object.test_quantity = "1 km"
    assert isinstance(test_object.test_quantity, pint.Quantity)
    assert test_object.test_quantity.check("[length]")
    assert test_object.test_quantity == ureg("1000 m")

    test_object.test_quantity = ureg("1 cm")
    assert isinstance(test_object.test_quantity, pint.Quantity)
    assert test_object.test_quantity.check("[length]")
    assert test_object.test_quantity.m == 1
    assert test_object.test_quantity == ureg("0.01 m")

    with pytest.raises(ValueError):
        test_object.test_quantity = "1 second"

    with pytest.raises(ValueError):
        test_object.test_quantity = 1 * ureg.J

    with pytest.raises(ValueError):
        test_object.test_quantity = ureg("1 kg")

    with pytest.raises(ValueError):
        test_object.test_quantity = ureg.Quantity(1, "liter")


def test_disk_cache(monkeypatch, tmp_path):
    monkeypatch.setattr(util, "cache_dir", tmp_path)

    @disk_cache
    def func(x, refresh_cache=False):
        return x + 1

    assert not os.path.isfile(os.path.join(tmp_path, "func(1).cache"))
    assert func(1) == 2
    assert os.path.isfile(os.path.join(tmp_path, "func(1).cache"))
    assert func(1) == 2
    with open(os.path.join(tmp_path, "func(1).cache"), "w+") as fp:
        json.dump(100, fp)
    assert func(1) == 100


def test_fsolve():
    assert fsolve(lambda x: x - 5, 0) == 5
    assert fsolve(lambda x: x - 5, 1) == 5
    assert fsolve(lambda x: (x - 5) * (x + 5), 6) == 5
    assert fsolve(lambda x: (x - 5) * (x + 5), -3) == -5
    assert fsolve(lambda x: (x - 5) * (x + 5), 4, 6) == 5
    assert fsolve(lambda x: (x - 5) * (x + 5), -4, -6) == -5
