import csv
import io
import re
import urllib.request
from collections import OrderedDict
from fractions import Fraction
from .util import sanitize_energy
import numpy as np
from math import pi as π


try:
    from .. import _ureg, _HAS_PINT
except ImportError:
    _HAS_PINT = False
    _ureg = None


class StateRegistry(tuple):
    def __getitem__(self, key):
        if isinstance(key, int):
            return super().__getitem__(key)
        if isinstance(key, slice):
            return StateRegistry(super().__getitem__(key))
        elif isinstance(key, str):
            return next(state for state in self if key in state.term)
        else:
            raise TypeError('key must be integer, slice, or term string')

    def __repr__(self):
        repr = '{:d} States (\n'.format(len(self))
        if self.__len__() <= 6:
            for state in self:
                repr += (str(state) + '\n')
        else:
            for state in self[:3]:
                repr += (str(state) + '\n')
            repr += '...\n'
            for state in self[-3:]:
                repr += (str(state) + '\n')
        repr = repr[:-1] + ')'
        return repr


class State(OrderedDict):
    def __init__(self, USE_UNITS=False, ureg=None, **state):
        self.USE_UNITS = USE_UNITS and _HAS_PINT
        if ureg:
            self._ureg = ureg
        else:
            self._ureg = _ureg

        if not self.USE_UNITS:
            self._ureg['hbar'] = 1
            self._ureg['ε_0'] = 1/(4*π)
            self._ureg['c'] = 137.03599908356244

        if 'energy' in state:
            energy = state['energy']
        elif 'Level (Ry)' in state:
            if self.USE_UNITS:
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
        return 'State({:} {:}: {:.4g})'.format(self.valence, self.term, self.energy)

    @property
    def energy(self):
        return self['energy']

    @property
    def configuration(self):
        return self['configuration']

    @property
    def valence(self):
        match = re.match(r'^\S+\.(?P<valence>\S+)', self.configuration)
        if match is None:
            return ''
        else:
            return match['valence']

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
            return 'LS'
        if 'J1' in self:
            return 'JJ'
        if 'K' in self:
            return 'LK'
        return None

    @property
    def term(self):
        if self.coupling is None:
            return 'Ionization Limit'

        P = '*' if self.parity == -1 else ''
        if self.coupling == 'LS':
            return '{:g}{:}{:}{:}'.format(2*self.S + 1, L_inv[self.L], Fraction(self.J), P)
        if self.coupling == 'JJ':
            return '({:},{:}){:}{:}'.format(Fraction(self.J1), Fraction(self.J2), Fraction(self.J), P)
        if self.coupling == 'LK':
            return '{:g}[{:}]{:}{:}'.format(2*self.S + 1, Fraction(self.K), Fraction(self.J), P)

    @property
    def transitions_down(self):
        try:
            return self['transitions_down']
        except KeyError:
            return []

    @property
    def transitions_up(self):
        try:
            return self['transitions_up']
        except KeyError:
            return []

    @property
    def transitions(self):
        return self.transitions_down + self.transitions_up

    @property
    def down(self):
        return self.transitions_down

    @property
    def up(self):
        return self.transitions_up

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

    def scalar_polarizability(self, ω=0):
        if not hasattr(ω, '__len__'):
            ω = [ω]

        hbar = self._ureg['hbar']
        ε_0 = self._ureg['ε_0']
        c = self._ureg['c']

        ω0 = [(t.Ef - t.Ei)/hbar for t in self.up]
        Γ = [t.Gamma for t in self.up]
        deg = np.array([(2*t.f.J + 1)/(2*t.i.J + 1)
                        for t in self.up])
        if self.USE_UNITS:
            if len(ω0):
                ω0 = self._ureg.Quantity(
                    np.array([ω.magnitude for ω in ω0]), ω0[0].units)
            else:
                ω0 = self._ureg.Quantity([], _ureg['Eh/hbar'])

            if len(Γ):
                Γ = self._ureg.Quantity(
                    np.array([γ.magnitude for γ in Γ]), Γ[0].units)
            else:
                Γ = self._ureg.Quantity([], _ureg['Eh/hbar'])
        else:
            ω0 = np.array(ω0)
            Γ = np.array(Γ)

        RME2 = 3*π*ε_0*c**3 * ω0**-3 * deg * Γ

        RME2 = np.tile(RME2, (np.size(ω), 1))
        ω0_m = np.tile(ω0, (np.size(ω), 1))
        ω_m = np.tile(ω, (np.size(ω0), 1)).T
        Δ = ω0_m/(ω0_m**2 - ω_m**2)
        α0_up = np.sum((2/3)*RME2*Δ, axis=1)

        ω0 = [(t.Ef - t.Ei)/hbar for t in self.down]
        Γ = [t.Gamma for t in self.down]
        if self.USE_UNITS:
            if len(ω0):
                ω0 = self._ureg.Quantity(
                    np.array([ω.magnitude for ω in ω0]), ω0[0].units)
            else:
                ω0 = self._ureg.Quantity([], _ureg['Eh/hbar'])

            if len(Γ):
                Γ = self._ureg.Quantity(
                    np.array([γ.magnitude for γ in Γ]), Γ[0].units)
            else:
                Γ = self._ureg.Quantity([], _ureg['Eh/hbar'])
        else:
            ω0 = np.array(ω0)
            Γ = np.array(Γ)

        RME2 = 3*π*ε_0*c**3 * ω0**-3 * Γ

        RME2 = np.tile(RME2, (np.size(ω), 1))
        ω0_m = np.tile(ω0, (np.size(ω), 1))
        ω_m = np.tile(ω, (np.size(ω0), 1)).T
        Δ = ω0_m/(ω0_m**2 - ω_m**2)
        α0_down = np.sum((2/3)*RME2*Δ, axis=1)

        α0 = α0_up - α0_down
        return α0


def download_nist_states(atom):
    url = 'https://physics.nist.gov/cgi-bin/ASD/energy1.pl'
    values = {
        'spectrum': atom,
        'units': 2,  # energy units {0: cm^-1, 1: eV, 2: Ry}
        'format': 3,  # format {0: HTML, 1: ASCII, 2: CSV, 3: TSV}
        'multiplet_ordered': 1,  # energy ordred
        'term_out': 'on',  # output the term symbol string
        'conf_out': 'on',  # output the configutation string
        'level_out': 'on',  # output the energy level
        'unc_out': 0,  # uncertainty on energy
        'j_out': 'on',  # output the J level
        'g_out': 'on',  # output the g-factor
        'lande_out': 'off'  # output experimentally measured g-factor
    }

    get_postfix = urllib.parse.urlencode(values)
    with urllib.request.urlopen(url + '?' + get_postfix) as response:
        response = response.read()

    data = csv.DictReader(io.StringIO(response.decode()), dialect='excel-tab')

    states = [row for row in data]
    for state in states:
        del state[None]

    return states


def get_states(atom, USE_UNITS=False, ureg=None):
    data = download_nist_states(atom)
    states = StateRegistry(
        State(**row, USE_UNITS=USE_UNITS, ureg=ureg) for row in data)
    return states


LS_term = re.compile(r'^(?P<S>\d+)(?P<L>[A-Z])\*?$')
JJ_term = re.compile(r'^\((?P<J1>\d+/?\d*),(?P<J2>\d+/?\d*)\)\*?$')
LK_term = re.compile(r'^(?P<S>\d+)\[(?P<K>\d+/?\d*)\]\*?$')

L = {
    'S': 0, 'P': 1, 'D': 2, 'F': 3,
    'G': 4, 'H': 5, 'I': 6, 'K': 7,
    'L': 8, 'M': 9, 'N': 10, 'O': 11,
    'Q': 12, 'R': 13, 'T': 14, 'U': 15}
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
