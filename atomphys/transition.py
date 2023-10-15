from collections import UserList
from collections.abc import Iterable
from math import inf
from math import pi as π
from typing import Any

import pint

from . import _ureg, state
from .util import default_units, fsolve


class Transition:

    _ureg: pint.UnitRegistry
    __matrix_element: pint.Quantity
    __type: str
    __atom = None
    __state_i = None
    __state_f = None

    def __init__(self, state_i=None, state_f=None, ureg=None, atom=None, **transition):
        self.__atom = atom

        self._ureg = _ureg
        if ureg is not None:
            self._ureg = ureg
        if atom is not None:
            self._ureg = atom._ureg

        if not state_i:
            state_i = state.State()
        if isinstance(state_i, dict):
            states = self.__atom.states.match(**state_i)
            if "energy" in state_i:
                self.__state_i = states(state_i["energy"])
            elif "En" in state_i:
                self.__state_i = states(state_i["En"])
            else:
                self.__state_i = states[0]
        else:
            self.__state_i = state_i
        self.__state_i.transitions.append(self)

        if not state_f:
            state_f = state.State()
        if isinstance(state_f, dict):
            states = self.__atom.states.match(**state_f)
            if "energy" in state_f:
                self.__state_f = states(state_f["energy"])
            elif "En" in state_f:
                self.__state_f = states(state_f["En"])
            else:
                self.__state_f = states[0]
        else:
            self.__state_f = state_f
        self.__state_f.transitions.append(self)

        self._matrix_element = 0
        for attr in [
            "energy",
            "En",
            "angular_frequency",
            "omega",
            "ω",
            "frequency",
            "nu",
            "ν",
            "wavelength",
            "λ",
            "matrix_element",
            "d",
            "Gamma",
            "Γ",
            "A",
        ]:
            if attr in transition:
                setattr(self, "_" + attr, transition[attr])

        if "type" in transition:
            self.__type = transition["type"]
        else:
            self.__type = ""

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
    def _energy(self):
        return self.f.energy - self.i.energy

    @_energy.setter
    @default_units("E_h")
    def _energy(self, energy):
        self.f._energy = self.i.energy + energy

    @property
    def energy(self):
        return self._energy

    @property
    def _En(self):
        return self._energy

    @_En.setter
    def _En(self, energy):
        self._energy = energy

    @property
    def En(self):
        return self._En

    @property
    def _angular_frequency(self):
        return (self._energy) / self._ureg.ħ

    @_angular_frequency.setter
    @default_units("THz")
    def _angular_frequency(self, omega):
        self._energy = omega * self._ureg.ħ

    @property
    def angular_frequency(self):
        return self._angular_frequency

    @property
    def _omega(self):
        return self._angular_frequency

    @_omega.setter
    def _omega(self, omega):
        self._angular_frequency = omega

    @property
    def omega(self):
        return self._omega

    @property
    def _ω(self):
        return self._angular_frequency

    @_ω.setter
    def _ω(self, omega):
        self._angular_frequency = omega

    @property
    def ω(self):
        return self._ω

    @property
    def _frequency(self):
        return (self._energy) / self._ureg.planck_constant

    @_frequency.setter
    @default_units("THz")
    def _frequency(self, frequency):
        self._energy = frequency * self._ureg.planck_constant

    @property
    def frequency(self):
        return self._frequency

    @property
    def _nu(self):
        return self._frequency

    @_nu.setter
    def _nu(self, frequency):
        self._frequency = frequency

    @property
    def nu(self):
        return self._nu

    @property
    def _ν(self):
        return self._frequency

    @_ν.setter
    def _ν(self, frequency):
        self._frequency = frequency

    @property
    def ν(self):
        return self._ν

    @property
    def _wavelength(self):
        try:
            return (self._ureg.planck_constant * self._ureg.c) / self._energy
        except ZeroDivisionError:
            return self._ureg.Quantity(inf, "nm")

    @_wavelength.setter
    @default_units("nm")
    def _wavelength(self, wavelength):
        self._energy = (self._ureg.planck_constant * self._ureg.c) / wavelength

    @property
    def wavelength(self):
        return self._wavelength

    @property
    def _λ(self):
        return self._wavelength

    @_λ.setter
    def _λ(self, wavelength):
        self._wavelength = wavelength

    @property
    def λ(self):
        return self._λ

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
    def _d(self):
        return self.__matrix_element

    @_d.setter
    def _d(self, matrix_element):
        self._matrix_element = matrix_element

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
    def _Γ(self):
        return self._Gamma

    @_Γ.setter
    def _Γ(self, Gamma):
        self._Gamma = Gamma

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
    def branching_ratio(self):
        return (self.Γ * self.f.τ).m_as("dimensionless")

    @property
    def saturation_intensity(self):
        h = self._ureg.planck_constant
        c = self._ureg.c
        return π * h * c * self.Γ / (3 * self.λ ** 3)

    @property
    def Isat(self):
        return self.saturation_intensity

    @property
    def σ0(self):
        ħ = self._ureg.ħ
        return ħ * self.ω * self.Γ / (2 * self.Isat)

    @property
    def cross_section(self):
        return self.σ0

    def polarizability(self, mJ_i=None, mJ_f=None, **kwargs):
        return self.__state_f.polarizability(
            mJ=mJ_f, **kwargs
        ) - self.__state_i.polarizability(mJ=mJ_i, **kwargs)

    @property
    def α(self):
        return self.polarizability

    @default_units("nm")
    def magic_wavelength(self, estimate, mJ_i=None, mJ_f=None, **kwargs):
        return fsolve(
            lambda λ: self.polarizability(λ=λ, mJ_i=None, mJ_f=None, **kwargs), estimate
        )

    @property
    def λ_magic(self):
        return self.magic_wavelength


class TransitionRegistry(UserList):
    _ureg: pint.UnitRegistry = None

    def __init__(self, data=[], ureg=None, atom=None):
        if not all(isinstance(transition, Transition) for transition in data):
            raise TypeError("TransitionRegistry can only contain transitions")
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

            try:
                quantity = self._ureg.Quantity(key)
                return self(quantity)
            except (pint.errors.UndefinedUnitError, pint.errors.DimensionalityError):
                pass

            raise KeyError(f"no transition {key} found")
        elif isinstance(key, float):
            wavelength = self._ureg.Quantity(key, "nm")
            return min(self, key=lambda t: abs(t.wavelength - wavelength))
        elif isinstance(key, self._ureg.Quantity):
            if key.check("[length]"):
                return min(self, key=lambda t: abs(t.wavelength - key))
            elif key.check("1/[time]"):
                return min(self, key=lambda t: abs(t.frequency - key))
            elif key.check("[energy]"):
                return min(self, key=lambda t: abs(t.energy - key))
        elif isinstance(key, Iterable):
            return TransitionRegistry([self(item) for item in key], ureg=self._ureg)
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

    def _assert_Transition(self, transition: Any):
        if not isinstance(transition, Transition):
            raise TypeError("TransitionRegistry can only contain transitions")

    def __setitem__(self, index: int, transition: Transition):
        self._assert_Transition(transition)
        super().__setitem__(index, transition)

    def insert(self, index: int, transition: Transition):
        self._assert_Transition(transition)
        super().insert(index, transition)

    def append(self, transition: Transition):
        self._assert_Transition(transition)
        super().append(transition)

    def extend(self, transitions):
        [self._assert_Transition(transition) for transition in transitions]
        super().extend(transitions)

    def up_from(self, state):
        return TransitionRegistry(
            [transition for transition in self if transition.i == state],
            ureg=self._ureg,
        )

    def down_from(self, state):
        return TransitionRegistry(
            [transition for transition in self if transition.f == state],
            ureg=self._ureg,
        )

    def to_dict(self):
        return [transition.to_dict() for transition in self]
