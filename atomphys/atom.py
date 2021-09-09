from .states import State, StateRegistry
from .transitions import Transition, TransitionRegistry
from .data import fetch_states, fetch_transitions, nuclear
from .util import parse_atom_name, parse_nuc_data
from math import pi as π
import os
import json
import math
import re

try:
    from . import _ureg, _HAS_PINT
except ImportError:
    _HAS_PINT = False
    _ureg = None

current_file = os.path.realpath(__file__)
directory = os.path.dirname(current_file)
periodic_table = os.path.join(directory, "data", "PeriodicTableJSON.json")
with open(periodic_table) as f:
    pt = json.load(f)
    symbols = [element['symbol'] for element in pt['elements']]

nuc_periodic_table = os.path.join(directory, "data", "NuclearPeriodicTableJSON.json")
with open(nuc_periodic_table) as f:
    nuc_ptable = json.load(f)


class Atom:

    name = ''
    nuc_spin = 0
    g_nuc = 0
    abundance = 1.0
    half_life = ''
    mass = 0.0

    def __init__(self, atom, USE_UNITS=True, ureg=None):
        self._USE_UNITS = USE_UNITS and _HAS_PINT

        if ureg:
            self._ureg = ureg
        else:
            self._ureg = _ureg

        try:
            self.load(atom)
        except FileNotFoundError:
            (
                atom_sym,
                atom_charge,
                atom_isotope,
            ) = parse_atom_name(atom)
            self.name = atom_sym + atom_charge
            self.isotope = atom_isotope
            self.load_nist(self.name)
            self.load_nuc(self.isotope + atom_sym)

        # reverse sort by Gamma first
        self._transitions.sort(key=lambda transition: transition.Gamma, reverse=True)
        # then sort by upper state energy
        self._transitions.sort(key=lambda transition: transition.Ef)
        # sort then by lower state energy
        self._transitions.sort(key=lambda transition: transition.Ei)
        # because sort is stable, this produces a list sorted by both upper and
        # lower state energy

        # index the transitions to the states
        for transition in self._transitions:
            transition._atom = self
            transition._state_i = next(state for state in self._states if state.energy == transition.Ei)
            transition._state_f = next(state for state in self._states if state.energy == transition.Ef)

        # index the states to the transitions
        for state in self._states:
            state._transitions_down = self._transitions.down_from(state)
            state._transitions_up = self._transitions.up_from(state)

    def __getitem__(self, state):
        return self.states[state]

    def __call__(self, state):
        return self.states[state]

    def __repr__(self):
        return (
            f'Isotope: {self.isotope}, I = {self.nuc_spin}\n'
            f'Ground State: {self.states[0]}\n'
            f'{len(self.states)} States\n'
            f'{len(self.transitions)} Transitions'
        )

    def to_dict(self):
        return {
            'name': self.name,
            'isotope': self.isotope,
            'I': self.I,
            'gI': self.gI,
            'abundance': self.abundance,
            'half_life': self.half_life,
            'mass': self.half_life,
            'states': self.states.to_dict(),
            'transitions': self.transitions.to_dict(),
        }

    def save(self, filename):
        with open(filename, 'w') as file:
            json.dump(self.to_dict(), file, indent=4, ensure_ascii=False)

    def load(self, filename):
        with open(filename) as file:
            data = json.load(file)

        self.name = data['name']
        self.isotope = data['isotope']
        self.nuc_spin = data['I']
        self.g_nuc = data['gI']
        self.abundance = data['abundance']
        self.half_life = data['half_life']
        self.mass = data['mass']

        self._states = StateRegistry(
            (State(**state, USE_UNITS=self._USE_UNITS, ureg=self._ureg) for state in data['states']), parent=self
        )
        self._transitions = TransitionRegistry(
            Transition(**transition, USE_UNITS=self._USE_UNITS, ureg=self._ureg) for transition in data['transitions']
        )
        print('loaded %s' % filename)

    def load_nist(self, name):
        if name in symbols:
            atom = name + ' i'
        elif name[-1] == '+' and name[:-1] in symbols:
            atom = name[:-1] + ' ii'
        else:
            atom = name
            raise Exception(f'{atom} does not match a known neutral atom or ionic ion name')

        self._states = StateRegistry(
            (State(**state, USE_UNITS=self._USE_UNITS, ureg=self._ureg, parent=self) for state in fetch_states(atom)),
            parent=self,
        )
        self._transitions = TransitionRegistry(
            Transition(**transition, USE_UNITS=self._USE_UNITS, ureg=self._ureg)
            for transition in fetch_transitions(atom)
        )

    def load_nuc(self, name):
        nuc_data = [row for row in nuc_ptable if name in row['Nuclide']][0]
        nuc_data = parse_nuc_data(nuc_data)
        self.nuc_spin = nuc_data['spin']
        self.g_nuc = nuc_data['gI']
        self.abundance = nuc_data['abundance']
        self.half_life = nuc_data['half_life']
        self.mass = nuc_data['mass']

    @property
    def states(self):
        return self._states

    @property
    def transitions(self):
        return self._transitions

    @property
    def units(self):
        return self._ureg

    @property
    def I(self):
        return self.nuc_spin

    @property
    def gI(self):
        return self.g_nuc
