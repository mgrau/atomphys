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
    __states: StateRegistry = []
    __transitions: TransitionRegistry = []

    def __init__(self, atom=None, ureg=None, refresh_cache=False):
        self._ureg = ureg if ureg is not None else _ureg

        if atom is None:
            return

        try:
            self.load(atom)
        except FileNotFoundError:
            self.load_nist(atom, refresh_cache)

        self.__transitions.sort()
        for state in self.__states:
            state.transitions.sort()

    def __call__(self, key):
        try:
            return self.__states(key)
        except KeyError:
            pass
        try:
            return self.__transitions(key)
        except KeyError:
            pass
        raise KeyError(f"no property of atom {self.name} {key}")

    def __repr__(self):
        ground_state = (
            f"Ground State: {self.__states[0].term}\n" if self.__states else ""
        )
        return (
            f"{ground_state}{len(self.__states)} States\n"
            f"{len(self.__transitions)} Transitions"
        )

    def _load_states(self, states):
        self.__states = StateRegistry(
            sorted([State(**state, atom=self) for state in states]),
            atom=self,
        )

    def _load_transitions(self, transitions):
        self.__transitions = TransitionRegistry(
            [Transition(**transition, atom=self) for transition in transitions],
            atom=self,
        )

    def to_dict(self):
        return {
            "name": self.name,
            "states": self.__states.to_dict(),
            "transitions": self.__transitions.to_dict(),
        }

    def save(self, filename):
        with open(filename, "w") as file:
            json.dump(self.to_dict(), file, indent=4, ensure_ascii=False)

    def load(self, filename):
        with open(filename) as file:
            data = json.load(file)

        self.name = data["name"]
        self._load_states(data["states"])
        self._load_transitions(data["transitions"])

    def load_nist(self, name, refresh_cache=False):
        if name in symbols:
            atom = name + " i"
        elif name[-1] == "+" and name[:-1] in symbols:
            atom = name[:-1] + " ii"
        elif name[-1] == "+" and name[-2].isdigit() and name[:-2] in symbols:
            atom = name[:-2] + " " + "i" * (1 + int(name[-2]))
        else:
            atom = name
            raise ValueError(
                f"{atom} does not match a known neutral atom or ionic ion name"
            )
        atom = atom.lower()

        self.name = name
        self._load_states(nist.parse_states(nist.fetch_states(atom, refresh_cache)))
        self._load_transitions(
            nist.parse_transitions(nist.fetch_transitions(atom, refresh_cache))
        )

    @property
    def states(self) -> StateRegistry:
        """StateRegistry: the atomic states."""
        return self.__states

    @property
    def transitions(self) -> TransitionRegistry:
        """TransitionRegistry: the atomic transitions."""
        return self.__transitions

    @property
    def units(self):
        """pint.UnitRegistry(): readonly access to the pint UnitRegistry used by the atom."""
        return self._ureg
