from collections import UserList
from collections.abc import Iterable
from math import inf
from math import pi as π

import pint

from . import _ureg
from .util import default_units, fsolve


class TransitionRegistry(UserList):
    __atom = None

    def __init__(self, data=[], atom=None):
        super().__init__(data)
        self.__atom = atom

    def __call__(self, key):
        if isinstance(key, int):
            return self[key]
        elif isinstance(key, str):
            try:
                quantity = self.__atom._ureg.Quantity(key)
                return self(quantity)
            except (pint.errors.UndefinedUnitError, pint.errors.DimensionalityError):
                pass

            try:
                return next(
                    transition
                    for transition in self
                    if (transition.i.match(key) or transition.f.match(key))
                )
            except StopIteration:
                pass

            try:
                for splitter in [
                    ":",
                    "to",
                    ",",
                    "<--->",
                    "<-->",
                    "<->",
                    "--->",
                    "-->",
                    "->",
                ]:
                    if splitter in key:
                        state_i, state_f = key.split(splitter)
                        return next(
                            transition
                            for transition in self
                            if (
                                (
                                    transition.i.match(state_i.strip())
                                    and transition.f.match(state_f.strip())
                                )
                                or (
                                    transition.i.match(state_f.strip())
                                    and transition.f.match(state_i.strip())
                                )
                            )
                        )
            except StopIteration:
                pass

            raise KeyError(f"no transition {key} found")
        elif isinstance(key, float):
            wavelength = self.__atom._ureg.Quantity(key, "nm")
            return min(self, key=lambda t: abs(t.wavelength - wavelength))
        elif isinstance(key, self.__atom._ureg.Quantity):
            if key.check("[length]"):
                return min(self, key=lambda t: abs(t.wavelength - key))
            elif key.check("1/[time]"):
                return min(self, key=lambda t: abs(t.frequency - key))
            elif key.check("[energy]"):
                return min(self, key=lambda t: abs(t.energy - key))
        elif isinstance(key, Iterable):
            return TransitionRegistry((self(item) for item in key), atom=self.__atom)
        else:
            raise TypeError("key must be integer index, term string")

    def __repr__(self):
        repr = f"{len(self)} Transitions (\n"
        if self.__len__() <= 6:
            for transition in self:
                repr += str(transition) + "\n"
        else:
            for transition in self[:3]:
                repr += str(transition) + "\n"
            repr += "...\n"
            for transition in self[-3:]:
                repr += str(transition) + "\n"
        repr = repr[:-1] + ")"
        return repr

    def up_from(self, state):
        return TransitionRegistry(
            [transition for transition in self if transition.i == state],
            atom=self.__atom,
        )

    def down_from(self, state):
        return TransitionRegistry(
            [transition for transition in self if transition.f == state],
            atom=self.__atom,
        )

    def to_dict(self):
        return [transition.to_dict() for transition in self]


class Transition:

    _ureg: pint.UnitRegistry
    __matrix_element: pint.Quantity
    __type: None
    __atom = None
    __state_i = None
    __state_f = None

    def __init__(self, state_i, state_f, ureg=None, atom=None, **transition):
        self.__atom = atom

        self._ureg = _ureg
        if ureg is not None:
            self._ureg = ureg
        if atom is not None:
            self._ureg = atom._ureg

        if isinstance(state_i, dict):
            states = self.__atom.states.match(**state_i)
            if "energy" in state_i:
                self.__state_i = states(state_i["energy"])
            elif "En" in state_i:
                self.__state_i = states(state_i["En"])
            else:
                self.__state_i = next(states)
        else:
            self.__state_i = state_i
        self.__state_i.transitions.append(self)

        if isinstance(state_f, dict):
            states = self.__atom.states.match(**state_f)
            if "energy" in state_f:
                self.__state_f = states(state_f["energy"])
            elif "En" in state_f:
                self.__state_f = states(state_f["En"])
            else:
                self.__state_f = next(states)
        else:
            self.__state_f = state_f
        self.__state_f.transitions.append(self)

        if "matrix_element" in transition:
            self._matrix_element = transition["matrix_element"]
        if "d" in transition:
            self._matrix_element = transition["d"]

        if "Gamma" in transition:
            self._Gamma = transition["Gamma"]
        if "Γ" in transition:
            self._Gamma = transition["Γ"]
        if "A" in transition:
            self._A = transition["A"]

        if "type" in transition:
            self.__type = transition["type"]

    def __repr__(self):
        return (
            f"Transition({self.i.name} <--> {self.f.name} : "
            f"λ={self.λ.to('nm'):0.5g~P}, "
            f"Γ=2π×{(self.Γ/(2*π)).to('Hz').to_compact():0.3g~P})"
        )

    def __lt__(self, other):
        return self.energy < other.energy

    def to_dict(self):
        return {
            "state_i": {"energy": f"{self.i.energy: ~P}", "term": self.i.term},
            "state_f": {"energy": f"{self.f.energy: ~P}", "term": self.f.term},
            "wavelength": f'{self.λ.to("nm"):0.3f~P}',
            "matrix_element": f'{self.d.to("e a0"):~P}',
            "type": self.type,
        }

    # ------
    # States
    # ------

    @property
    def state_i(self):
        return self.__state_i

    @property
    def i(self):
        return self.__state_i

    @property
    def state_f(self):
        return self.__state_f

    @property
    def f(self):
        return self.__state_f

    # ------
    # Energy
    # ------

    @property
    def energy(self):
        return self.f.energy - self.i.energy

    @property
    def angular_frequency(self):
        ħ = self._ureg.ħ
        return (self.energy) / ħ

    @property
    def ω(self):
        return self.angular_frequency

    @property
    def frequency(self):
        return self.ω / (2 * π)

    @property
    def ν(self):
        return self.frequency

    @property
    def wavelength(self):
        c = self._ureg.c
        try:
            return c / self.ν
        except ZeroDivisionError:
            return inf

    @property
    def λ(self):
        return self.wavelength

    # --------------
    # matrix element
    # --------------

    @property
    def _matrix_element(self):
        return self.__matrix_element

    @_matrix_element.setter
    @default_units("e a0")
    def _matrix_element(self, matrix_element):
        self.__matrix_element = matrix_element

    @property
    def _Gamma(self):
        ε_0 = self._ureg.ε_0
        ħ = self._ureg.ħ
        c = self._ureg.c

        ω = self.ω
        J = self.f.J
        d = self.matrix_element
        return (ω ** 3) / ((3 * π * ε_0 * ħ * c ** 3) * (2 * J + 1)) * d ** 2

    @_Gamma.setter
    @default_units("s^-1")
    def _Gamma(self, Γ):
        ε_0 = self._ureg.ε_0
        ħ = self._ureg.ħ
        c = self._ureg.c

        ω = self.ω
        J = self.f.J
        self._matrix_element = (
            (3 * π * ε_0 * ħ * c ** 3) / (ω ** 3) * (2 * J + 1) * Γ
        ) ** (1 / 2)

    @property
    def _A(self):
        return self._Gamma

    @_A.setter
    def _A(self, A):
        self._Gamma = A

    @property
    def matrix_element(self):
        return self.__matrix_element

    @property
    def d(self):
        return self.__matrix_element

    @property
    def Gamma(self):
        return self._Gamma

    @property
    def Γ(self):
        return self._Gamma

    @property
    def A(self):
        return self._Gamma

    @property
    def type(self):
        return self.__type

    # ---------------------
    # high level properties
    # ---------------------

    @property
    def saturation_intensity(self):
        h = self._ureg.h
        c = self._ureg.c
        return π * h * c * self.Γ / (3 * self.λ ** 3)

    @property
    def Isat(self):
        return self.saturation_intensity

    @property
    def branching_ratio(self):
        return (self.Γ * self.f.τ).m_as("dimensionless")

    @property
    def σ0(self):
        ħ = self._ureg.ħ
        return ħ * self.ω * self.Γ / (2 * self.Isat)

    @property
    def cross_section(self):
        return self.σ0

    def magic_wavelength(self, estimate, mJ_i=None, mJ_f=None, **kwargs):
        α_i = self._state_i.α
        α_f = self._state_f.α

        def f(λ):
            return α_i(mJ=mJ_i, λ=λ, **kwargs) - α_f(mJ=mJ_f, λ=λ, **kwargs)

        return fsolve(f, estimate)

    @property
    def λ_magic(self):
        return self.magic_wavelength
