import csv
import io
import re
import enum
import urllib.request
from fractions import Fraction
from collections.abc import Iterable
from .util import sanitize_energy
from .data import nist
from .calc import polarizability
from .transitions import TransitionRegistry
from .laser import Laser
from .constants import gs
from math import pi as π
from itertools import chain


try:
    from . import _ureg, _HAS_PINT
except ImportError:
    _HAS_PINT = False
    _ureg = None


class Coupling(enum.Enum):
    LS = 'LS'  # Russell-Saunders
    jj = 'jj'
    LK = 'LK'  # pair coupling


class StateRegistry(list):
    _parent = None

    def __init__(self, *args, parent=None, **kwargs):
        super().__init__(*args, **kwargs)
        self._parent = parent

    def __getitem__(self, key):
        if isinstance(key, int):
            return super().__getitem__(key)
        elif isinstance(key, slice):
            return StateRegistry(super().__getitem__(key), parent=self._parent)
        elif isinstance(key, str):
            return next(state for state in self if state.match(key))
        elif isinstance(key, Iterable):
            return StateRegistry((self.__getitem__(item) for item in key), parent=self._parent)
        elif isinstance(key, float):
            energy = self._parent._ureg.Quantity(
                key, 'E_h') if self._parent.USE_UNITS else key
            return min(self, key=lambda state: abs(state.energy - energy))
        elif self._parent.USE_UNITS and isinstance(key, self._parent._ureg.Quantity):
            return min(self, key=lambda state: abs(state.energy - key))
        else:
            raise TypeError('key must be integer, slice, or term string')

    def __call__(self, key):
        return self.__getitem__(key)

    def __repr__(self):
        repr = f'{len(self)} States (\n'
        if self.__len__() <= 6:
            for state in self:
                repr += f'{state}\n'
        else:
            for state in self[:3]:
                repr += f'{state}\n'
            repr += '...\n'
            for state in self[-3:]:
                repr += f'{state}\n'
        repr = repr[:-1] + ')'
        return repr

    def __add__(self, other):
        assert isinstance(other, StateRegistry)
        return StateRegistry(list(self) + list(other), parent=self._parent)

    def to_dict(self):
        return [state.to_dict() for state in self]


class State(dict):

    _USE_UNITS = False
    _ureg = {}
    _transitions_down = []
    _transitions_up = []

    def __init__(self, USE_UNITS=False, ureg=None, **state):
        self._USE_UNITS = USE_UNITS and _HAS_PINT
        if ureg and self._USE_UNITS:
            self._ureg = ureg
        elif self._USE_UNITS:
            self._ureg = _ureg
        else:
            self._ureg = {}

        if not self._USE_UNITS:
            self._ureg['ħ'] = 1
            self._ureg['ε_0'] = 1/(4*π)
            self._ureg['c'] = 137.03599908356244

        if 'energy' in state:
            if self._USE_UNITS:
                energy = self._ureg.Quantity(state['energy'])
            else:
                energy = float(state['energy'])
        elif 'Level (Ry)' in state:
            if self._USE_UNITS:
                energy = self._ureg.Quantity(
                    float(sanitize_energy(state['Level (Ry)'])), 'Ry').to('Eh')
            else:
                energy = 0.5 * float(sanitize_energy(state['Level (Ry)']))
        else:
            energy = 0.0

        if 'configuration' in state:
            configuration = state['configuration']
        elif 'Configuration' in state:
            configuration = state['Configuration']
        else:
            configuration = ''

        if 'J' in state:
            try:
                J = float(Fraction(state['J'].strip('?')))
            except:
                J = 0
        else:
            J = 0

        if 'term' in state:
            term = parse_term(state['term'])
        elif 'Term' in state:
            term = parse_term(state['Term'])
        else:
            term = {}

        super(State, self).__init__(
            {'energy': energy, 'configuration': configuration, 'J': J,  **term})

    def __repr__(self):
        fmt = '0.4g~P' if self._USE_UNITS else '0.4g'
        return f'State({self.name}: {self.energy:{fmt}})'

    def to_dict(self):
        return {'energy': str(self.energy), 'configuration': self.configuration, 'term': self.term, 'J': str(Fraction(self.J))}

    def match(self, name):
        return name in self.name

    @property
    def energy(self):
        return self['energy']

    @property
    def configuration(self):
        return self['configuration']

    @property
    def valence(self):
        match = re.match(r'^(\S*\.)?(?P<valence>\S+)', self.configuration)
        if match is None:
            return ''
        else:
            return match['valence']

    @property
    def name(self):
        return f'{self.valence} {self.term}'

    @property
    def J(self):
        return self['J']

    @property
    def S(self):
        return self['S']

    @property
    def L(self):
        return self['L']

    @property
    def J1(self):
        return self['J1']

    @property
    def J2(self):
        return self['J2']

    @property
    def K(self):
        return self['K']

    @property
    def parity(self):
        return self['parity']

    @property
    def coupling(self):
        if 'L' in self:
            return Coupling.LS
        if 'J1' in self:
            return Coupling.jj
        if 'K' in self:
            return Coupling.LK
        return None

    @property
    def g(self):
        if self.coupling == Coupling.LS:
            L, S, J = self.L, self.S, self.J
            return (gs+1)/2 + (gs-1)/2 * (S*(S+1) - L*(L-1))/(J*(J+1))
        else:
            return None

    @property
    def term(self):
        if self.coupling is None:
            return 'Ionization Limit'

        P = '*' if self.parity == -1 else ''
        if self.coupling == Coupling.LS:
            return f'{2*self.S + 1:g}{L_inv[self.L]}{Fraction(self.J)}{P}'
        if self.coupling == Coupling.jj:
            return f'({Fraction(self.J1)},{Fraction(self.J2)}){Fraction(self.J)}{P}'
        if self.coupling == Coupling.LK:
            return f'{2*self.S + 1:g}[{Fraction(self.K)}]{Fraction(self.J)}{P}'

    @property
    def transitions_down(self):
        try:
            return self._transitions_down
        except KeyError:
            return []

    @property
    def transitions_up(self):
        try:
            return self._transitions_up
        except KeyError:
            return []

    @property
    def transitions(self):
        return TransitionRegistry(self.transitions_down + self.transitions_up, parent=self)

    @property
    def to(self):
        return self.transitions

    @property
    def down(self):
        return TransitionRegistry(self.transitions_down, parent=self)

    @property
    def up(self):
        return TransitionRegistry(self.transitions_up, parent=self)

    @property
    def lifetime(self):
        Gamma = [transition.Gamma for transition in self.transitions_down]
        try:
            lifetime = 1/sum(Gamma)
        except ZeroDivisionError:
            lifetime = float('inf')

        return lifetime

    @property
    def τ(self):
        return self.lifetime

    def polarizability(self, mJ=None, laser=None, **kwargs):
        if laser is None:
            laser = Laser(**kwargs)
        else:
            laser = Laser(laser=laser, **kwargs)
        return polarizability.total(self, mJ=mJ,
                                    omega=laser.omega,
                                    A=laser.A, theta_k=laser.theta_k, theta_p=laser.theta_p)

    @property
    def α(self):
        return self.polarizability


LS_term = re.compile(r'^(?P<S>\d+)(?P<L>[A-Z])\*?')
JJ_term = re.compile(r'^\((?P<J1>\d+/?\d*),(?P<J2>\d+/?\d*)\)\*?')
LK_term = re.compile(r'^(?P<S>\d+)\[(?P<K>\d+/?\d*)\]\*?')

L = {
    'S': 0, 'P': 1, 'D': 2, 'F': 3,
    'G': 4, 'H': 5, 'I': 6, 'K': 7,
    'L': 8, 'M': 9, 'N': 10, 'O': 11,
    'Q': 12, 'R': 13, 'T': 14, 'U': 15,
    'V': 16, 'W': 17, 'X': 18, 'Y': 19}
L_inv = {value: key for key, value in L.items()}


def parse_term(term):
    '''
    parse term symbol string
    '''
    if term == 'Limit':
        return {}

    parity = -1 if '*' in term else 1

    match = LS_term.match(term)
    if match is None:
        match = JJ_term.match(term)
    if match is None:
        match = LK_term.match(term)
    if match is None:
        return {'partiy': parity}

    def convert(key, value):
        if key == 'S':
            return float(Fraction((int(value)-1)/2))
        if key == 'J1' or key == 'J2' or key == 'K':
            return float(Fraction(value))
        if key == 'L':
            return L[value]

    term = {key: convert(key, value)
            for key, value in match.groupdict().items()}

    return {**term, 'parity': parity}
