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


class State:

    _ureg: pint.UnitRegistry
    __energy: pint.Quantity
    __quantum_numbers: dict
    __atom = None
    __transitions: TransitionRegistry

    def __init__(self, term=None, energy=None, ureg=None, atom=None, **state):
        self.__atom = atom

        self._ureg = _ureg
        if ureg is not None:
            self._ureg = ureg
        if atom is not None:
            self._ureg = atom._ureg

        if energy:
            self._energy = energy
        elif "En" in state:
            self._En = state["En"]
        else:
            self._energy = 0

        self.__quantum_numbers = {}
        if term:
            self.__quantum_numbers = parse_term(term)

        for qN in ["S", "L", "J1", "J2", "S2", "K", "J", "n", "parity"]:
            if qN in state:
                self.__quantum_numbers[qN] = state[qN]
        if "J" not in self.__quantum_numbers:
            self.J = 0

        self.__transitions = TransitionRegistry(atom=self.__atom)

    def __repr__(self):
        name = f"{self.name}: " if self.name else ""
        return f"State({name}{self.energy:0.4g~P})"

    def __getattr__(self, name: str) -> Any:
        if name in self.__quantum_numbers:
            return self.__quantum_numbers[name]
        raise AttributeError(f"{self.__class__.__name__} has no attribute {name}")

    def __eq__(self, other):
        return (
            self.energy == other.energy
            and self.quantum_numbers == other.quantum_numbers
        )

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
    def term(self):
        return print_term(**self.__quantum_numbers)

    @property
    def name(self):
        if self.term is None:
            return None
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
        try:
            L, S, J = self.L, self.S, self.J
        except AttributeError:
            return None
        if J == 0:
            return 0
        return (J * (J + 1) - S * (S + 1) + L * (L + 1)) / (2 * J * (J + 1)) + gs * (
            J * (J + 1) + S * (S + 1) - L * (L + 1)
        ) / (2 * J * (J + 1))

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

    def polarizability(self, laser=None, mJ=None, **kwargs):
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


class StateRegistry(UserList):
    _ureg: pint.UnitRegistry = None

    def __init__(self, data=[], ureg=None, atom=None):
        if not all(isinstance(state, State) for state in data):
            raise TypeError("StateRegistry can only contain states")
        super().__init__(data)

        if atom:
            self._ureg = atom._ureg
        elif ureg:
            self._ureg = ureg
        else:
            self._ureg = _ureg

    def __call__(self, key):
        if isinstance(key, int):
            return self[key]
        elif isinstance(key, str):
            try:
                return self(self._ureg.Quantity(key))
            except (pint.errors.UndefinedUnitError, pint.errors.DimensionalityError):
                pass

            try:
                return next(state for state in self if state.match(key))
            except StopIteration:
                pass

            raise KeyError(f"no state {key} found")
        elif isinstance(key, float):
            energy = self._ureg.Quantity(key, "E_h")
            return min(self, key=lambda state: abs(state.energy - energy))
        elif isinstance(key, self._ureg.Quantity):
            return min(self, key=lambda state: abs(state.energy - key))
        elif isinstance(key, Iterable):
            return StateRegistry([self(item) for item in key], ureg=self._ureg)
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

    def _assert_State(self, state: Any):
        if not isinstance(state, State):
            raise TypeError("StateRegistry can only contain states")

    def __setitem__(self, index: int, state: State):
        self._assert_State(state)
        super().__setitem__(index, state)

    def insert(self, index: int, state: State):
        self._assert_State(state)
        super().insert(index, state)

    def append(self, state: State):
        self._assert_State(state)
        super().append(state)

    def extend(self, states):
        [self._assert_State(state) for state in states]
        super().extend(states)

    def search(self, func):
        def search_func(state):
            try:
                return func(state)
            except Exception:
                return False

        return StateRegistry(list(filter(search_func, self)), ureg=self._ureg)

    def match(self, **kwargs):
        kwargs.pop("energy", None)
        kwargs.pop("En", None)
        return self.search(
            lambda state: all(getattr(state, key) == val for key, val in kwargs.items())
        )

    def to_dict(self):
        return [state.to_dict() for state in self]
