import enum
from collections import UserList
from collections.abc import Iterable
from fractions import Fraction
from typing import Any

import pint

from . import _ureg
from .calc import polarizability
from .constants import gs
from .laser import Laser
from .term import L_inv, parse_term, print_term
from .transition import TransitionRegistry
from .util import default_units


class Coupling(enum.Enum):
    LS = "LS"  # Russell-Saunders
    jj = "jj"
    LK = "LK"  # pair coupling


class StateRegistry(UserList):
    __atom = None

    def __init__(self, data=[], atom=None):
        super().__init__(data)
        self.__atom = atom

    def __call__(self, key):
        if isinstance(key, int):
            return self[key]
        elif isinstance(key, str):
            try:
                return self(self.__atom._ureg.Quantity(key))
            except (pint.errors.UndefinedUnitError, pint.errors.DimensionalityError):
                pass

            try:
                return next(state for state in self if state.match(key))
            except StopIteration:
                pass

            raise KeyError(f"no state {key} found")
        elif isinstance(key, float):
            energy = self.__atom._ureg.Quantity(key, "E_h")
            return min(self, key=lambda state: abs(state.energy - energy))
        elif isinstance(key, self.__atom._ureg.Quantity):
            return min(self, key=lambda state: abs(state.energy - key))
        elif isinstance(key, Iterable):
            return StateRegistry([self(item) for item in key], atom=self.__atom)
        else:
            raise TypeError(
                "key must be integer index, term string, energy, or iterable"
            )

    def __repr__(self):
        repr = f"{len(self)} States (\n"
        if self.__len__() <= 6:
            for state in self:
                repr += f"{state}\n"
        else:
            for state in self[:3]:
                repr += f"{state}\n"
            repr += "...\n"
            for state in self[-3:]:
                repr += f"{state}\n"
        repr = repr[:-1] + ")"
        return repr

    def search(self, func):
        def search_func(state):
            try:
                return func(state)
            except BaseException:
                return False

        return StateRegistry(list(filter(search_func, self)), atom=self.__atom)

    def match(self, **kwargs):
        kwargs.pop("energy", None)
        kwargs.pop("En", None)
        return self.search(
            lambda state: all(getattr(state, key) == val for key, val in kwargs.items())
        )

    def to_dict(self):
        return [state.to_dict() for state in self]


class State:

    _ureg: pint.UnitRegistry
    __energy: pint.Quantity
    __quantum_numbers: dict = {}
    __configuration: str
    __atom = None
    __transitions: TransitionRegistry
    _transitions_down = []
    _transitions_up = []

    def __init__(self, ureg=None, atom=None, **state):
        self.__atom = atom

        self._ureg = _ureg
        if ureg is not None:
            self._ureg = ureg
        if atom is not None:
            self._ureg = atom._ureg

        self._energy = 0
        if "energy" in state:
            self._energy = state["energy"]
        if "En" in state:
            self._energy = state["En"]

        if "term" in state:
            self.__quantum_numbers = parse_term(state["term"])

        for qN in ["S", "L", "J", "J1", "J2", "K", "n", "parity"]:
            if qN in state:
                self.__quantum_numbers[qN] = state[qN]

        if "configuration" in state:
            self.__configuration = state["configuration"]

        self.__transitions = TransitionRegistry([], atom=self.__atom)

    def __repr__(self):
        try:
            name = f"{self.name}: "
        except RuntimeError:
            name = ""

        energy = f"{self.energy:0.4g~P}"
        return f"State({name}{energy})"

    def __getattr__(self, name: str) -> Any:
        if name in self.__quantum_numbers:
            return self.__quantum_numbers[name]
        raise AttributeError(f"{self.__class__.__name__} has not attribute {name}")

    def __lt__(self, other):
        return self.energy < other.energy

    def to_dict(self):
        return {
            "name": self.name,
            "energy": f"{self.energy:~P}",
            "term": self.term,
            "quantum numbers": self.quantum_numbers,
        }

    def match(self, name):
        return name.lower() in self.name.lower()

    # ------
    # Energy
    # ------

    @property
    def _energy(self) -> pint.Quantity:
        return self.__energy

    @_energy.setter
    @default_units("E_h")
    def _energy(self, energy):
        self.__energy = energy

    @property
    def energy(self) -> pint.Quantity:
        return self.__energy

    @property
    def _En(self) -> pint.Quantity:
        return self.__energy

    @_En.setter
    @default_units("E_h")
    def _En(self, energy):
        self.__energy = energy

    @property
    def En(self) -> pint.Quantity:
        return self.__energy

    @property
    def quantum_numbers(self) -> dict:
        return self.__quantum_numbers

    @property
    def configuration(self) -> str:
        return self.__configuration

    @property
    def term(self):
        return print_term(**self.__quantum_numbers)

    @property
    def name(self):
        if "n" in self.__quantum_numbers:
            return f"{self.n}{L_inv[self.L]}{Fraction(self.J)}, {self.term}"
        return f"{self.term}"

    # -----------
    # Transitions
    # -----------

    @property
    def transitions(self):
        return self.__transitions

    @property
    def to(self):
        return self.__transitions

    @property
    def transitions_down(self):
        return self.__transitions.down_from(self)

    @property
    def transitions_up(self):
        return self.__transitions.up_from(self)

    @property
    def down(self):
        return self.transitions_down

    @property
    def up(self):
        return self.transitions_up

    # ---------------------
    # high level properties
    # ---------------------

    @property
    def coupling(self):
        if all(key in self.__quantum_numbers.keys() for key in ["S", "L"]):
            return Coupling.LS
        if all(key in self.__quantum_numbers.keys() for key in ["J1", "J2"]):
            return Coupling.jj
        if all(key in self.__quantum_numbers.keys() for key in ["S2", "K"]):
            return Coupling.LK
        return None

    @property
    def g(self):
        if self.coupling == Coupling.LS:
            L, S, J = self.L, self.S, self.J
            return (gs + 1) / 2 + (gs - 1) / 2 * (S * (S + 1) - L * (L - 1)) / (
                J * (J + 1)
            )
        else:
            return None

    @property
    def lifetime(self):
        Gamma = sum([transition.Gamma for transition in self.transitions_down])
        if Gamma == 0:
            return self._ureg.Quantity(float("inf"), "s")
        else:
            return 1 / Gamma

    @property
    def τ(self):
        return self.lifetime

    def polarizability(self, mJ=None, laser=None, **kwargs):
        if laser is None:
            laser = Laser(**kwargs)
        else:
            laser = Laser(laser=laser, **kwargs)
        return polarizability.total(
            self,
            mJ=mJ,
            omega=laser.omega,
            A=laser.A,
            theta_k=laser.theta_k,
            theta_p=laser.theta_p,
        )

    @property
    def α(self):
        return self.polarizability
