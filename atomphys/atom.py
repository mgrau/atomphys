from .states import State, StateRegistry
from .transitions import Transition, TransitionRegistry
from .isotopes import Isotope, IsotopeRegistry
from .data import fetch_states, fetch_transitions, fetch_isotopes
from math import pi as π
import os
import json
import math
import re
from fractions import Fraction

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
    pt_symbols = [element["symbol"] for element in pt["elements"]]
    pt_names = [element["name"] for element in pt["elements"]]


class Atom:
    symbol = ""

    def __init__(self, atom, USE_UNITS=True, ureg=None):
        self._USE_UNITS = USE_UNITS and _HAS_PINT

        if ureg:
            self._ureg = ureg
        else:
            self._ureg = _ureg

        try:
            self.load(atom)
        except FileNotFoundError:
            self.name = ""
            self.isotope = ""
            self.nuclear_spin = ""

            atom_regex = re.search("^(\d*)([A-Za-z]+)([0-9+]*)", atom)
            if atom_regex is not None:
                atom_isotope = atom_regex.group(1)
                atom_sym = atom_regex.group(2)
                atom_charge = atom_regex.group(3)
            else:
                raise ValueError("regex did not yield a valid atom name string")
            self.symbol = atom_sym
            self.name = pt_names[pt_symbols.index(atom_sym)]
            self.load_nist(atom_sym + atom_charge)
            self.load_nuc(self.name, atom_sym, atom_isotope)

        self._isotopes = self._isotopes

        # reverse sort by Gamma first
        self._transitions.sort(key=lambda transition: transition.Gamma, reverse=True)
        # then sort by upper state energy
        self._transitions.sort(key=lambda transition: transition.Ef)
        # sort then by lower state energy
        self._transitions.sort(key=lambda transition: transition.Ei)
        # because sort is stable, this produces a list sorted by both upper and lower state energy

        # index the transitions to the states
        for transition in self._transitions:
            transition._atom = self
            transition._state_i = next(
                state for state in self._states if state.energy == transition.Ei
            )
            transition._state_f = next(
                state for state in self._states if state.energy == transition.Ef
            )

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
            f"Ground State: {self.states[0]}\n"
            f"{len(self.states)} States\n"
            f"{len(self.transitions)} Transitions\n"
            f"{len(self.isotopes)} Isotopes"
        )

    def to_dict(self):
        return {
            "symbol": self.symbol,
            "name": self.name,
            "nuclear_spin": self.nuclear_spin,
            "states": self.states.to_dict(),
            "transitions": self.transitions.to_dict(),
            "isotopes": self.isotopes.to_dict(),
        }

    def save(self, filename):
        with open(filename, "w") as file:
            json.dump(self.to_dict(), file, indent=4, ensure_ascii=False)

    def load(self, filename):
        with open(filename) as file:
            data = json.load(file)

        self.name = data["name"]
        self.symbol = data["symbol"]
        self._states = StateRegistry(
            (
                State(**state, USE_UNITS=self._USE_UNITS, ureg=self._ureg)
                for state in data["states"]
            ),
            parent=self,
        )
        self._transitions = TransitionRegistry(
            Transition(**transition, USE_UNITS=self._USE_UNITS, ureg=self._ureg)
            for transition in data["transitions"]
        )

    def load_nist(self, symbol):
        if symbol in pt_symbols:
            atom = symbol + " i"
        elif symbol[-1] == "+" and symbol[:-1] in pt_symbols:
            atom = symbol[:-1] + " ii"
        else:
            atom = symbol
            raise Exception(
                f"{atom} does not match a known neutral atom or ionic ion name"
            )

        self._states = StateRegistry(
            (
                State(**state, USE_UNITS=self._USE_UNITS, ureg=self._ureg)
                for state in fetch_states(atom)
            ),
            parent=self,
        )
        self._transitions = TransitionRegistry(
            Transition(**transition, USE_UNITS=self._USE_UNITS, ureg=self._ureg)
            for transition in fetch_transitions(atom)
        )

    def load_nuc(self, name, sym, iso):
        self._isotopes = IsotopeRegistry(
            (
                Isotope(**isotope, USE_UNITS=self._USE_UNITS, ureg=self._ureg)
                for isotope in fetch_isotopes(name, sym)
            ),
            parent=self,
        )

    @property
    def states(self):
        return self._states

    @property
    def transitions(self):
        return self._transitions

    @property
    def isotopes(self):
        return self._isotopes

    @property
    def units(self):
        return self._ureg
