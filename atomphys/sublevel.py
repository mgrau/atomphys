import csv
import io
import re
import enum
import urllib.request
from fractions import Fraction
from collections.abc import Iterable
from .data import nist
from .calc import polarizability
from .transitions import TransitionRegistry
from .laser import Laser
from .constants import gs
from math import nan, pi as π, inf
from itertools import chain
from fractions import Fraction
from .util import parse_time

try:
    from . import _ureg, _HAS_PINT
except ImportError:
    _HAS_PINT = False
    _ureg = None


class SublevelRegistry(list):
    _parent = None

    def __init__(self, *args, parent=None, **kwargs):
        super().__init__(*args, **kwargs)
        self._parent = parent

    def __getitem__(self, key):
        if isinstance(key, int):
            return super().__getitem__(key)
        elif isinstance(key, slice):
            return SublevelRegistry(
                super().__getitem__(key), parent=self._parent)
        elif isinstance(key, str):
            return next(isotope for isotope in self if isotope.match(key))
        elif isinstance(key, Iterable):
            return SublevelRegistry(
                (self.__getitem__(item) for item in key), parent=self._parent
            )
        else:
            raise TypeError("key must be integer, slice, or term string")

    def __call__(self, key):
        return self.__getitem__(key)

    def __repr__(self):
        repr = f"{len(self)} Isotopes (\n"
        if self.__len__() <= 6:
            for isotope in self:
                repr += f"{isotope}\n"
        else:
            for isotope in self[:3]:
                repr += f"{isotope}\n"
            repr += "...\n"
            for isotope in self[-3:]:
                repr += f"{isotope}\n"
        repr = repr[:-1] + ")"
        return repr

    def __add__(self, other):
        assert isinstance(other, SublevelRegistry)
        return SublevelRegistry(list(self) + list(other), parent=self._parent)

    def to_dict(self):
        return [isotope.to_dict() for isotope in self]


class Sublevel(dict):
    _USE_UNITS = False
    _ureg = {}

    def __init__(self, USE_UNITS=False, ureg=None, **isotope):
        self._USE_UNITS = USE_UNITS and _HAS_PINT
        if ureg and self._USE_UNITS:
            self._ureg = ureg
        elif self._USE_UNITS:
            self._ureg = _ureg
        else:
            self._ureg = {}

        if not self._USE_UNITS:
            self._ureg["ħ"] = 1
            self._ureg["ε_0"] = 1 / (4 * π)
            self._ureg["c"] = 137.03599908356244

        if "Nuclide" in isotope:
            nuclide = re.search("^\\d+[m]*\\d*[A-Za-z]+",
                                isotope["Nuclide"]).group()
        else:
            nuclide = ""
        
        super(Sublevel, self).__init__(
            {
                "I": I,
                "J": J,
            }
        )

    def __repr__(self):
        return f"isotope({self.nuclide}: atomic_mass = {self.atomic_mass}, spin = {self.spin}, abundance = {self.abundance}, half_life: {self.half_life})"

    def to_dict(self):
        return {
            "I": self.I,
            "J": self.J,
        }

    def match(self, nuclide):
        return nuclide in self.nuclide

    @property
    def I(self):
        return self["I"]

    @property
    def L(self):
        return self["L"]

    @property
    def S(self):
        return self["S"]

    @property
    def J(self):
        return self["J"]

    @property
    def F(self):
        return self["F"]
