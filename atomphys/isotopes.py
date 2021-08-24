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
            nuclide = re.search('^\d+[m]*\d*[A-Za-z]+',isotope['Nuclide']).group()
        else:
            nuclide = ''

        if 'Atomic number' in isotope:
            if isotope['Atomic number'] is not None:
                if len(isotope['Atomic number']) > 0:
                    try:
                        atomic_number = int(re.search('\d+',isotope['Atomic number']).group())
                    except AttributeError:
                        atomic_number = isotope['Atomic number']
                else:
                    atomic_number = 0.0
            else:
                atomic_number = 0.0
        else:
            atomic_number = None

        if 'Neutron number' in isotope:
            if isotope['Neutron number'] is not None:
                if len(isotope['Neutron number']) > 0:
                    try:
                        neutron_number = int(re.search('\d+',isotope['Neutron number']).group())
                    except AttributeError:
                        neutron_number = isotope['Neutron number']
                else:
                    neutron_number = 0.0
            else:
                neutron_number = 0.0
        else:
            neutron_number = None

        if 'Atomic mass' in isotope:
            if isotope['Atomic mass'] is not None:
                if len(isotope['Atomic mass']) > 0:
                    try:
                        atomic_mass = float(re.search('\d+\.*\d*',isotope['Atomic mass']).group())
                    except AttributeError:
                        atomic_mass = isotope['Atomic mass']
                else:
                    atomic_mass = 0.0
            else:
                atomic_mass = 0.0
        else:
            atomic_mass = None

        if 'Half-life' in isotope:
            if 'Stable' in isotope['Half-life']:
                half_life = self._ureg.Quantity(inf,'s')
            elif len(isotope['Half-life'])==0:
                half_life = None
            else:
                #half_life = isotope['Half-life']
                #print(isotope['Half-life'])
                s = isotope['Half-life']
                s = re.search('^(\d+\.*\d*).*\\xa0(\w+)$',s)
                if s is not None:
                    #print(s)
                    value = float(s.group(1))
                    t_factor = s.group(2)

                    if t_factor == 'h':
                        t_factor = 'hours'
                    elif t_factor == 'y':
                        t_factor = 'years'
                    elif t_factor == 'μs':
                        t_factor = 'microsecond'

                    if self._USE_UNITS:
                        try:
                            half_life = self._ureg.Quantity(value,t_factor)
                        except Exception:
                            half_life = value * parse_time(t_factor)
                    else:
                        half_life = value * parse_time(t_factor)
                else:
                    half_life = 0.0
        else:
            half_life = None

        if 'Spin (physics)' in isotope:
            if len(isotope['Spin (physics)']) > 0:
                try:
                    spin = Fraction(re.search('\d+[/]*\d*',isotope['Spin (physics)']).group())
                except AttributeError:
                    spin = isotope['Spin (physics)']
            else:
                spin = None
        else:
            spin = None

        if 'Natural abundance' in isotope:
            if len(isotope['Natural abundance']) > 0:
                try:
                    abundance = float(re.search('^\d+\.\d+',isotope['Natural abundance']).group())
                except AttributeError:
                    abundance = isotope['Natural abundance']
            else:
                abundance = None
        else:
            abundance = None

        if 'mag_moment_μN' in isotope:
            if isotope['mag_moment_μN'] is not None:
                try:
                    value = float(re.search('\d+\.+\d+',isotope['mag_moment_μN']).group())
                except AttributeError:
                    value = isotope['mag_moment_μN']

                if self._USE_UNITS:
                    try:
                        magnetic_moment = self._ureg.Quantity(value,'nuclear_magneton')
                    except Exception:
                        magnetic_moment = value 
                else:
                    magnetic_moment = value 

            else:
                magnetic_moment = 0.0
        else:
            magnetic_moment = None

        super(Isotope, self).__init__(
            {'nuclide': nuclide, 'atomic_number': atomic_number, 'neutron_number': neutron_number, 
            'atomic_mass': atomic_mass, 'spin': spin, 'magnetic_moment': magnetic_moment, 
            'half_life': half_life, 'abundance': abundance})

    def __repr__(self):
        return f'isotope({self.nuclide}: atomic_mass = {self.atomic_mass}, spin = {self.spin}, abundance = {self.abundance}, half_life: {self.half_life})'

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

    @property
    def Z(self):
        return self['atomic_mass']

    @property
    def A(self):
        return self['atomic_number']

    @property
    def N(self):
        return self['neutron_number']

    @property
    def nuclide(self):
        return self['nuclide']

    @property
    def spin(self):
        return self['spin']

    @property
    def half_life(self):
        return self['half_life']

    @property
    def abundance(self):
        return self['abundance']

    @property
    def atomic_mass(self):
        return self['atomic_mass']

    @property
    def neutron_number(self):
        return self['neutron_number']

    @property
    def proton_number(self):
        return self['proton_number']

    @property
    def magnetic_moment(self):
        return self['magnetic_moment']

    @property
    def μ(self):
        return self['magnetic_moment']

    @property
    def gI(self):
        return self['magnetic_moment']