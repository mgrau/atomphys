from .data import State, Transition, get_states, get_transitions
from math import pi as π
import os
import json


try:
    from . import _ureg, _HAS_PINT
except ImportError:
    _HAS_PINT = False
    _ureg = None

current_file = os.path.realpath(__file__)
directory = os.path.dirname(current_file)
periodic_table = os.path.join(
    directory, "data", "PeriodicTableJSON.json")
with open(periodic_table) as f:
    pt = json.load(f)
    symbols = [element['symbol'] for element in pt['elements']]


class Atom():

    def __init__(self, atom, USE_UNITS=True, ureg=None):
        self.USE_UNITS = USE_UNITS and _HAS_PINT

        if ureg:
            self._ureg = ureg
        else:
            self._ureg = _ureg

        if atom in symbols:
            atom = atom + ' i'
        elif atom[-1] == '+' and atom[:-1] in symbols:
            atom = atom[:-1] + ' ii'
        else:
            raise Exception(
                atom + ' does not match a known neutral atom or ionic ion name')

        self._states = get_states(
            atom, USE_UNITS=self.USE_UNITS, ureg=self._ureg)
        self._transitions = get_transitions(
            atom, USE_UNITS=self.USE_UNITS, ureg=self._ureg)

        # sort by upper state energy first
        self._transitions.sort(key=lambda transition: transition.Ef)
        # sort then by lower state energy
        # because sort is stable, this produces a list sorted by both upper and lower state energy
        self._transitions.sort(key=lambda transition: transition.Ei)

        # index the transitions to the states
        for transition in self._transitions:
            transition._state_i = next(
                state for state in self._states if state.energy == transition.Ei)
            transition._state_f = next(
                state for state in self._states if state.energy == transition.Ef)

        # index the states to the transitions
        for state in self._states:
            state._transitions_down = self._transitions.down_from(state)
            state._transitions_up = self._transitions.up_from(state)

    def __getitem__(self, state):
        return self.states[state]

    def __repr__(self):
        repr = ''
        repr += 'Ground State: {:}\n'.format(self.states[0])
        repr += '{:d} States\n'.format(len(self.states))
        repr += '{:d} Transitions'.format(len(self.transitions))
        return repr

    @property
    def states(self):
        return self._states

    @property
    def transitions(self):
        return self._transitions
