import enum
from collections.abc import Iterable
from fractions import Fraction
from typing import Any

import pint

from . import _ureg
from .calc import polarizability
from .constants import gs
from .laser import Laser
from .term import parse_term, print_term
from .transition import TransitionRegistry
from .util import default_units


class Coupling(enum.Enum):
    LS = "LS"  # Russell-Saunders
    jj = "jj"
    LK = "LK"  # pair coupling


class StateRegistry(list):
    _parent = None

    def __init__(self, *args, parent=None, **kwargs):
        super().__init__(*args, **kwargs)
        self._parent = parent

    def __getitem__(self, key):
        if isinstance(key, int):
            return super().__getitem__(key)
        elif isinstance(key, slice):
            return StateRegistry(super().__getitem__(key), parent=self._parent)
        elif isinstance(key, str):
            return next(state for state in self if state.match(key))
        elif isinstance(key, Iterable):
            return StateRegistry(
                (self.__getitem__(item) for item in key), parent=self._parent
            )
        elif isinstance(key, float):
            energy = self._parent._ureg.Quantity(key, "E_h")
            return min(self, key=lambda state: abs(state.energy - energy))
        elif isinstance(key, self._parent._ureg.Quantity):
            return min(self, key=lambda state: abs(state.energy - key))
        else:
            raise TypeError("key must be integer, slice, or term string")

    def __call__(self, key):
        return self.__getitem__(key)

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

    def __add__(self, other):
        assert isinstance(other, StateRegistry)
        return StateRegistry(list(self) + list(other), parent=self._parent)

    def to_dict(self):
        return [state.to_dict() for state in self]


class State:

    _ureg: pint.UnitRegistry
    __energy: pint.Quantity
    __quantum_numbers: dict = {}
    __atom = None
    __transitions: TransitionRegistry
    _transitions_down = []
    _transitions_up = []

    def __init__(self, ureg=None, atom=None, **state):
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

        for qN in ["S", "L", "J", "J1", "J2", "K", "parity"]:
            if qN in state:
                self.__quantum_numbers[qN] = state["qN"]

        if "term" in state:
            self.__quantum_numbers = {
                **self.__quantum_numbers,
                **parse_term(state["term"]),
            }

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

    def to_dict(self):
        return {
            "energy": str(self.energy),
            "configuration": self.configuration,
            "term": self.term,
            "J": str(Fraction(self.J)),
        }

    def match(self, name):
        return name in self.name

    # ------
    # Energy
    # ------

    @property
    def _energy(self):
        return self.__energy

    @_energy.setter
    @default_units("Ry")
    def _energy(self, energy):
        self.__energy = energy

    @property
    def energy(self):
        return self.__energy

    @property
    def _En(self):
        return self.__energy

    @_En.setter
    @default_units("Ry")
    def _En(self, energy):
        self.__energy = energy

    @property
    def En(self):
        return self.__energy

    @property
    def quantum_numbers(self):
        return self.__quantum_numbers

    # @property
    # def valence(self):
    #     match = re.match(r"^(\S*\.)?(?P<valence>\S+)", self.configuration)
    #     if match is None:
    #         return ""
    #     else:
    #         return match["valence"]

    @property
    def name(self):
        return f"{self.term}"

    @property
    def term(self):
        return print_term(**self.__quantum_numbers)

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
    def transitions_down(self):
        try:
            return self._transitions_down
        except KeyError:
            return []

    @property
    def transitions_up(self):
        try:
            return self._transitions_up
        except KeyError:
            return []

    @property
    def transitions(self):
        return TransitionRegistry(
            self.transitions_down + self.transitions_up, parent=self
        )

    @property
    def to(self):
        return self.transitions

    @property
    def down(self):
        return TransitionRegistry(self.transitions_down, parent=self)

    @property
    def up(self):
        return TransitionRegistry(self.transitions_up, parent=self)

    @property
    def lifetime(self):
        Gamma = [transition.Gamma for transition in self.transitions_down]
        try:
            lifetime = 1 / sum(Gamma)
        except ZeroDivisionError:
            lifetime = float("inf") * self.__units.s

        return lifetime

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
