from .states import State, StateRegistry
from .transitions import Transition, TransitionRegistry
from .data import fetch_states, fetch_transitions, nuclear
from .util import parse_atom_name, parse_nuc_data
from math import pi as π
from math import isclose
import os
import json
import math
import re
import csv

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
            if len(atom_isotope) > 0:
                self.load_nuc(self.isotope + atom_sym)

        self._transitions._sort()
        self._transitions.index_to_states()
        self._states.index_to_transitions()

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
            'mass': self.mass,
            'states': self.states.to_dict(),
            'transitions': self.transitions.to_dict(),
        }

    def save(self, filename):
        with open(filename, 'w') as file:
            json.dump(self.to_dict(), file, indent=4, ensure_ascii=False)

    def save_csv(self):
        states_file = self.isotope + self.name + '_states.csv'
        transitions_file = self.isotope + self.name + '_transitions.csv'

        states_data = self.states.to_dict()
        col_headers = list(states_data[0].keys())
        with open(states_file, 'w') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=col_headers)
            writer.writeheader()
            [writer.writerow(data) for data in states_data]

        transitions_data = self.transitions.to_dict()
        col_headers = list(transitions_data[0].keys())
        with open(transitions_file, 'w') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=col_headers)
            writer.writeheader()
            [writer.writerow(data) for data in transitions_data]

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

    def load_csv(self):
        states_file = self.isotope + self.name + '_states.csv'
        transitions_file = self.isotope + self.name + '_transitions.csv'

        states_data = []
        with open(states_file, mode='r') as f:
            dict_reader = csv.DictReader(f)
            for row in dict_reader:
                states_data.append(row)

        transitions_data = []
        with open(transitions_file, mode='r') as f:
            dict_reader = csv.DictReader(f)
            for row in dict_reader:
                transitions_data.append(row)

        self._states = StateRegistry(
            (State(**state, USE_UNITS=self._USE_UNITS, ureg=self._ureg) for state in states_data), parent=self
        )
        self._transitions = TransitionRegistry(
            Transition(**transition, USE_UNITS=self._USE_UNITS, ureg=self._ureg) for transition in transitions_data
        )
        print('loaded states and transitions from csv')

    def load_nist(self, name):
        if name in symbols:
            atom = name + ' i'
        elif name[-1] == '+' and name[:-1] in symbols:
            atom = name[:-1] + ' ii'
        else:
            atom = name
            raise Exception(f'{atom} does not match a known neutral atom or ionic ion name')

        self._states = StateRegistry(
            (State(**state, USE_UNITS=self._USE_UNITS, ureg=self._ureg, atom=self) for state in fetch_states(atom)),
            atom=self,
        )
        self._transitions = TransitionRegistry(
            (
                Transition(**transition, USE_UNITS=self._USE_UNITS, ureg=self._ureg, atom=self)
                for transition in fetch_transitions(atom)
            ),
            atom=self,
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
