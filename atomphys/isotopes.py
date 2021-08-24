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
from math import pi as π, inf
from itertools import chain

try:
    from . import _ureg, _HAS_PINT
except ImportError:
    _HAS_PINT = False
    _ureg = None


class IsotopeRegistry(list):
    _parent = None

    def __init__(self, *args, parent=None, **kwargs):
        super().__init__(*args, **kwargs)
        self._parent = parent

    def __getitem__(self, key):
        if isinstance(key, int):
            return super().__getitem__(key)
        elif isinstance(key, slice):
            return IsotopeRegistry(super().__getitem__(key), parent=self._parent)
        elif isinstance(key, str):
            return next(isotope for isotope in self if isotope.match(key))
        elif isinstance(key, Iterable):
            return IsotopeRegistry((self.__getitem__(item) for item in key), parent=self._parent)
        else:
            raise TypeError('key must be integer, slice, or term string')

    def __call__(self, key):
        return self.__getitem__(key)

    def __repr__(self):
        repr = f'{len(self)} Isotopes (\n'
        if self.__len__() <= 6:
            for isotope in self:
                repr += f'{isotope}\n'
        else:
            for isotope in self[:3]:
                repr += f'{isotope}\n'
            repr += '...\n'
            for isotope in self[-3:]:
                repr += f'{isotope}\n'
        repr = repr[:-1] + ')'
        return repr

    def __add__(self, other):
        assert isinstance(other, IsotopeRegistry)
        return IsotopeRegistry(list(self) + list(other), parent=self._parent)

    def to_dict(self):
        return [isotope.to_dict() for isotope in self]


class Isotope(dict):
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
            self._ureg['ħ'] = 1
            self._ureg['ε_0'] = 1 / (4 * π)
            self._ureg['c'] = 137.03599908356244

        if 'Nuclide' in isotope:
            nuclide = isotope['Nuclide']
        else:
            nuclide = ''

        if 'Atomic number' in isotope:
            atomic_number = isotope['Atomic number']
        else:
            atomic_number = None

        if 'Neutron number' in isotope:
            neutron_number = isotope['Neutron number']
        else:
            neutron_number = None

        if 'Atomic mass' in isotope:
            atomic_mass = isotope['Atomic mass']
        else:
            atomic_mass = None

        if 'Half-life' in isotope:
            half_life = isotope['Half-life']
        else:
            half_life = inf

        if 'Spin (physics)' in isotope:
            spin = isotope['Spin (physics)']
        else:
            spin = 0.0

        if 'Natural abundance' in isotope:
            abundance = isotope['Natural abundance']
        else:
            abundance = 1.0

        super(Isotope, self).__init__(
            {'nuclide': nuclide, 'spin': spin, 'half-life': half_life, 'abundance': abundance})

    def __repr__(self):
        fmt = '0.4g~P' if self._USE_UNITS else '0.4g'
        return f'isotope({self.nuclide}: {self.spin:{fmt}})'

    def to_dict(self):
        return {'nuclide': str(self.nuclide), 'atomic_number': self.atomic_number,
                'neutron_number': self.neutron_number, 'atomic_mass': self.atomic_mass,
                'half_life': self.half_life, 'abundance': self.abundance,
                'spin': str(Fraction(self.spin))}

    def match(self, nuclide):
        return nuclide in self.nuclide

    @property
    def I(self):
        return self['spin']
