import csv
import io
import urllib.request
from collections.abc import Iterable
from .util import sanitize_energy, fsolve, dipole_allowed
from .data import nist
from itertools import combinations

from math import pi as π
from math import inf


try:
    from . import _ureg, _HAS_PINT
except ImportError:
    _HAS_PINT = False
    _ureg = None


class TransitionRegistry(list):
    _atom = None

    def __init__(self, *args, atom=None, **kwargs):
        super().__init__(*args, **kwargs)
        self._atom = atom

    def __getitem__(self, key):
        if isinstance(key, int):
            return super().__getitem__(key)
        elif isinstance(key, slice):
            return TransitionRegistry(super().__getitem__(key), atom=self._atom)
        elif isinstance(key, str):
            if ':' in key:
                state_i, state_f = key.split(':')
                return next(
                    transition for transition in self if ((transition.i.match(state_i) and transition.f.match(state_f)))
                )
            return next(
                transition
                for transition in self
                if (
                    (transition.i.match(key) and transition.i != self._atom)
                    or (transition.f.match(key) and transition.f != self._atom)
                )
            )
        elif isinstance(key, Iterable):
            return TransitionRegistry((self.__getitem__(item) for item in key), atom=self._atom)
        elif isinstance(key, float):
            energy = self._atom._ureg.Quantity(key, self._atom._energy_unit) if self._atom._USE_UNITS else key
            return min(
                self, key=lambda transition: min(abs(transition.i.energy - energy), abs(transition.f.energy - energy))
            )
        elif self._atom._USE_UNITS and isinstance(key, self._atom._ureg.Quantity):
            return min(self, key=lambda transition: min(abs(transition.i.energy - key), abs(transition.f.energy - key)))
        else:
            raise TypeError('key must be integer, slice, or term string')

    def __call__(self, key):
        return self.__getitem__(key)

    def __repr__(self):
        repr = f'{len(self)} Transitions (\n'
        if self.__len__() <= 6:
            for transition in self:
                repr += str(transition) + '\n'
        else:
            for transition in self[:3]:
                repr += str(transition) + '\n'
            repr += '...\n'
            for transition in self[-3:]:
                repr += str(transition) + '\n'
        repr = repr[:-1] + ')'
        return repr

    def __add__(self, other):
        assert isinstance(other, TransitionRegistry)
        return TransitionRegistry(list(self) + list(other), atom=self._atom)

    def up_from(self, state):
        return TransitionRegistry((transition for transition in self if transition.i == state), atom=self._atom)

    def down_from(self, state):
        return TransitionRegistry((transition for transition in self if transition.f == state), atom=self._atom)

    def to_dict(self):
        return [transition.to_dict() for transition in self]

    def populate_allowed(self, E_lower_max=0.12, wl_min=350, wl_max=2000):
        energy_max = (_ureg['h'] * _ureg['c'] / _ureg.Quantity(wl_min,'nm')).to(self._atom._energy_unit)
        energy_min = (_ureg['h'] * _ureg['c'] / _ureg.Quantity(wl_max,'nm')).to(self._atom._energy_unit)

        for pair in combinations(self._atom.states, 2):
            pair = sorted(list(pair),key=lambda state: state.energy)
            si,sf, = pair[0:2]
            if (
                dipole_allowed(si, sf)
                and si.energy <= _ureg.Quantity(E_lower_max, 'E_h')
                and not any([(si, sf) == (t.i, t.f) for t in self])
                and energy_min < abs(sf.energy - si.energy) < energy_max
            ):
                self.append(Transition(Ei=si.energy, Ef=sf.energy, USE_UNITS=True, ureg=_ureg, atom=self._atom))
        self._sort()
        self.index_to_states()
        self._atom._states.index_to_transitions()
                

    def _sort(self):
        # reverse sort by Gamma first
        self.sort(key=lambda transition: transition.Gamma, reverse=True)
        # then sort by upper state energy
        self.sort(key=lambda transition: transition.Ef)
        # sort then by lower state energy
        self.sort(key=lambda transition: transition.Ei)
        # because sort is stable, this produces a list sorted by both upper and
        # lower state energy

class Transition(dict):

    _USE_UNITS = False
    _ureg = {}
    _atom = None
    _state_i = None
    _state_f = None

    def __init__(self, USE_UNITS=False, ureg=None, atom=None, **transition):
        self._USE_UNITS = USE_UNITS and _HAS_PINT
        if ureg and self._USE_UNITS:
            self._ureg = ureg
        elif self._USE_UNITS:
            self._ureg = _ureg
        else:
            self._ureg = {}

        if not self._USE_UNITS:
            self._ureg['ħ'] = 1
            self._ureg['h'] = 2 * π
            self._ureg['ε_0'] = 1 / (4 * π)
            self._ureg['c'] = 137.03599908356244
        
        self._atom = atom

        if 'Gamma' in transition:
            if self._USE_UNITS:
                Gamma = self._ureg.Quantity(transition['Gamma']).to(self._atom._linewidth_unit)
            else:
                Gamma = float(transition['Gamma'])
        elif 'Aki(s^-1)' in transition:
            if self._USE_UNITS:
                Gamma = self._ureg.Quantity(float(transition['Aki(s^-1)']), 's^-1').to(self._atom._linewidth_unit)
            else:
                Gamma = 2.4188843265856806e-17 * float(transition['Aki(s^-1)'])
        else:
            if self._USE_UNITS:
                Gamma = self._ureg.Quantity(0.0, self._atom._linewidth_unit)
            else:
                Gamma = 0.0

        if 'Ei' in transition:
            if self._USE_UNITS:
                Ei = self._ureg.Quantity(transition['Ei']).to(self._atom._energy_unit)
            else:
                Ei = float(transition['Ei'])
        elif 'Ei(Ry)' in transition:
            if self._USE_UNITS:
                Ei = self._ureg.Quantity(float(sanitize_energy(transition['Ei(Ry)'])), 'Ry').to(self._atom._energy_unit)
            else:
                Ei = 0.5 * float(sanitize_energy(transition['Ei(Ry)']))
        elif 'Ei(cm-1)' in transition:
            if self._USE_UNITS:
                Ei = self._ureg.Quantity(float(sanitize_energy(transition['Ei(cm-1)'])), 'h*c/cm').to(self._atom._energy_unit)
            else:
                Ei = 0.5 * float(sanitize_energy(transition['Ei(cm-1)']))
        else:
            Ei = 0.0

        if 'Ef' in transition:
            if self._USE_UNITS:
                Ef = self._ureg.Quantity(transition['Ef']).to(self._atom._energy_unit)
            else:
                Ef = float(transition['Ef'])
        elif 'Ek(Ry)' in transition:
            if self._USE_UNITS:
                Ef = self._ureg.Quantity(float(sanitize_energy(transition['Ek(Ry)'])), 'Ry').to(self._atom._energy_unit)
            else:
                Ef = 0.5 * float(sanitize_energy(transition['Ek(Ry)']))
        elif 'Ek(cm-1)' in transition:
            if self._USE_UNITS:
                Ef = self._ureg.Quantity(float(sanitize_energy(transition['Ek(cm-1)'])), 'h*c/cm').to(self._atom._energy_unit)
            else:
                Ef = 0.5 * float(sanitize_energy(transition['Ek(cm-1)']))
        else:
            Ef = 0.0

        super(Transition, self).__init__({'Ei': Ei, 'Ef': Ef, 'Gamma': Gamma})

    def __repr__(self):
        fmt = '0.10g~P' if self._USE_UNITS else '0.10g'
        if self.i is not None:
            state_i = f'{self.i.valence} {self.i.term}'
        else:
            state_i = f'{self.Ei:{fmt}}'
        if self.f is not None:
            state_f = f'{self.f.valence} {self.f.term}'
        else:
            state_f = f'{self.Ef:{fmt}}'

        return f'Transition({state_i} <---> {state_f}, ' f'λ={self.λ_nm:{fmt}}, ' f'Γ=2π×{self.Γ_MHz/(2*π):{fmt}})'

    def to_dict(self):
        fmt = '0.10g~P' if self._USE_UNITS else '0.10g'
        return {
            'i': self._atom.states.index(self.i),
            'i_config': self.i.configuration,
            'i_term': self.i.term,
            'f': self._atom.states.index(self.f),
            'f_config': self.f.configuration,
            'f_term': self.f.term,
            'Ei': f'{self.Ei.to(self._atom._energy_unit):{fmt}}',
            'Ef': f'{self.Ef.to(self._atom._energy_unit):{fmt}}',
            'Gamma': f'{self.Gamma.to(self._atom._linewidth_unit):{fmt}}',
        }

    def index_to_states(self):
        self._state_i = next(state for state in self._atom._states if state.energy == self.Ei)
        self._state_f = next(state for state in self._atom._states if state.energy == self.Ef)

    @property
    def Ei(self):
        return self['Ei']

    @property
    def Ef(self):
        return self['Ef']

    @property
    def Gamma(self):
        return self['Gamma']

    @property
    def Γ(self):
        return self['Gamma']

    @property
    def Gamma_MHz(self):
        return self.Γ_MHz

    @property
    def Γ_MHz(self):
        if self._USE_UNITS:
            return self.Γ.to('MHz')
        else:
            return self.Γ * 41341373335.18245  # E_h / ħ / MHz

    @property
    def i(self):
        return self._state_i

    @property
    def f(self):
        return self._state_f

    @property
    def ω(self):
        ħ = self._ureg['ħ']
        return (self.Ef - self.Ei) / ħ

    @property
    def angular_frequency(self):
        return self.ω

    @property
    def ν(self):
        return self.ω / (2 * π)

    @property
    def frequency(self):
        return self.ν

    @property
    def λ(self):
        c = self._ureg['c']
        try:
            return c / self.ν
        except ZeroDivisionError:
            return inf

    @property
    def wavelength(self):
        return self.λ

    @property
    def λ_nm(self):
        if self._USE_UNITS:
            return self.λ.to('nm')
        else:
            return self.λ * 0.052917721090397746  # nm / a_0

    @property
    def wavelength_nm(self):
        return self.λ_nm

    @property
    def saturation_intensity(self):
        h = self._ureg['h']
        c = self._ureg['c']
        return π * h * c * self.Γ / (3 * self.λ ** 3)

    @property
    def Isat(self):
        return self.saturation_intensity

    @property
    def branching_ratio(self):
        r = self.Γ * self.f.τ
        if self._USE_UNITS and isinstance(r, self._ureg.Quantity):
            r = r.m
        return r

    @property
    def reduced_dipole_matrix_element(self):
        ω0 = self.ω
        ε_0 = self._ureg['ε_0']
        ħ = self._ureg['ħ']
        c = self._ureg['c']
        Jg = self.i.J
        Je = self.f.J
        Γ = self.Γ
        return ((3 * π * ε_0 * ħ * c ** 3) / (ω0 ** 3) * (2 * Je + 1) * Γ) ** (1 / 2)

    @property
    def reduced_dipole_matrix_element_conjugate(self):
        Jg = self.i.J
        Je = self.f.J
        d = self.reduced_dipole_matrix_element
        return (-1) ** (Je - Jg) * d

    @property
    def d(self):
        return self.reduced_dipole_matrix_element

    @property
    def d_conj(self):
        return self.reduced_dipole_matrix_element_conjugate

    @property
    def σ0(self):
        ħ = self._ureg['ħ']
        return ħ * self.ω * self.Γ / (2 * self.Isat)

    @property
    def cross_section(self):
        return self.σ0

    def magic_wavelength(self, estimate, mJ_i=None, mJ_f=None, **kwargs):
        c = self._ureg['c']
        α_i = self._state_i.α
        α_f = self._state_f.α

        def f(λ):
            return α_i(mJ=mJ_i, λ=λ, **kwargs) - α_f(mJ=mJ_f, λ=λ, **kwargs)

        return fsolve(f, estimate)

    @property
    def λ_magic(self):
        return self.magic_wavelength
