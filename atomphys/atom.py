import json
import os

import pint

from . import _ureg
from .data import nist
from .state import State, StateRegistry
from .transition import Transition, TransitionRegistry

current_file = os.path.realpath(__file__)
directory = os.path.dirname(current_file)
periodic_table = os.path.join(directory, "data", "PeriodicTableJSON.json")
with open(periodic_table) as f:
    pt = json.load(f)
    symbols = [element["symbol"] for element in pt["elements"]]


class Atom:
    """An atom object, containing states and transitions

    Attributes:
        name (str): The name of the atom
    """

    _ureg: pint.UnitRegistry
    name: str = ""

    def __init__(self, atom, ureg=None):
        self._ureg = ureg if ureg is not None else _ureg

        try:
            self.load(atom)
        except FileNotFoundError:
            self.name = atom
            self.load_nist(self.name)

        # reverse sort by Gamma first
        # self._transitions.sort(key=lambda transition: transition.Gamma, reverse=True)
        # then sort by upper state energy
        # self._transitions.sort(key=lambda transition: transition.Ef)
        # sort then by lower state energy
        # self._transitions.sort(key=lambda transition: transition.Ei)
        # because sort is stable, this produces a list sorted by both upper
        # and lower state energy

        # index the transitions to the states
        # for transition in self._transitions:
        #     transition._atom = self
        #     transition._state_i = next(
        #         state for state in self._states if state.energy == transition.Ei
        #     )
        #     transition._state_f = next(
        #         state for state in self._states if state.energy == transition.Ef
        #     )

        # index the states to the transitions
        # for state in self._states:
        #     state._transitions_down = self._transitions.down_from(state)
        #     state._transitions_up = self._transitions.up_from(state)

    def __getitem__(self, state):
        return self.states[state]

    def __call__(self, state):
        return self.states[state]

    def __repr__(self):
        return (
            f"Ground State: {self.states[0]}\n"
            f"{len(self.states)} States\n"
            f"{len(self.transitions)} Transitions"
        )

    def to_dict(self):
        return {
            "name": self.name,
            "states": self.states.to_dict(),
            "transitions": self.transitions.to_dict(),
        }

    def save(self, filename):
        with open(filename, "w") as file:
            json.dump(self.to_dict(), file, indent=4, ensure_ascii=False)

    def load(self, filename):
        with open(filename) as file:
            data = json.load(file)

        self.name = data["name"]
        self._states = StateRegistry(
            (State(**state, ureg=self._ureg) for state in data["states"]),
            parent=self,
        )
        self._transitions = TransitionRegistry(
            Transition(**transition, ureg=self._ureg)
            for transition in data["transitions"]
        )

    def load_nist(self, name):
        if name in symbols:
            atom = name + " i"
        elif name[-1] == "+" and name[:-1] in symbols:
            atom = name[:-1] + " ii"
        else:
            atom = name
            raise Exception(
                f"{atom} does not match a known neutral atom or ionic ion name"
            )

        self._states = StateRegistry(
            (
                State(**state, atom=self)
                for state in nist.parse_states(nist.fetch_states(atom))
            ),
            parent=self,
        )
        # self._transitions = TransitionRegistry(
        #     Transition(**transition, ureg=self._ureg)
        #     for transition in fetch_transitions(atom)
        # )

    @property
    def states(self) -> StateRegistry:
        """StateRegistry: the atomic states."""
        return self._states

    @property
    def transitions(self) -> TransitionRegistry:
        """TransitionRegistry: the atomic transitions."""
        return self._transitions

    @property
    def units(self):
        """pint.UnitRegistry(): readonly access to the pint UnitRegistry used by the atom."""
        return self._ureg
