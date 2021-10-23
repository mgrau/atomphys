from collections.abc import Iterable
from math import inf
from math import pi as π

import pint

from . import _ureg
from .data.nist import remove_annotations
from .util import fsolve


class TransitionRegistry(list):
    _parent = None

    def __init__(self, *args, parent=None, **kwargs):
        super().__init__(*args, **kwargs)
        self._parent = parent

    def __getitem__(self, key):
        if isinstance(key, int):
            return super().__getitem__(key)
        elif isinstance(key, slice):
            return TransitionRegistry(super().__getitem__(key), parent=self._parent)
        elif isinstance(key, str):
            if ":" in key:
                state_i, state_f = key.split(":")
                return next(
                    transition
                    for transition in self
                    if ((transition.i.match(state_i) and transition.f.match(state_f)))
                )
            return next(
                transition
                for transition in self
                if (
                    (transition.i.match(key) and transition.i != self._parent)
                    or (transition.f.match(key) and transition.f != self._parent)
                )
            )
        elif isinstance(key, Iterable):
            return TransitionRegistry(
                (self.__getitem__(item) for item in key), parent=self._parent
            )
        elif isinstance(key, float):
            energy = self._parent._ureg.Quantity(key, "E_h")
            return min(
                self,
                key=lambda transition: min(
                    abs(transition.i.energy - energy),
                    abs(transition.f.energy - energy),
                ),
            )
        elif self._parent._USE_UNITS and isinstance(key, self._parent._ureg.Quantity):
            return min(
                self,
                key=lambda transition: min(
                    abs(transition.i.energy - key),
                    abs(transition.f.energy - key),
                ),
            )
        else:
            raise TypeError("key must be integer, slice, or term string")

    def __call__(self, key):
        return self.__getitem__(key)

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

    def __add__(self, other):
        assert isinstance(other, TransitionRegistry)
        return TransitionRegistry(list(self) + list(other), parent=self._parent)

    def up_from(self, state):
        return TransitionRegistry(
            (transition for transition in self if transition.i == state),
            parent=self._parent,
        )

    def down_from(self, state):
        return TransitionRegistry(
            (transition for transition in self if transition.f == state),
            parent=self._parent,
        )

    def to_dict(self):
        return [transition.to_dict() for transition in self]


class Transition(dict):

    _ureg: pint.UnitRegistry
    __matrix_element: pint.Quantity
    __type: None
    __atom = None
    __state_i = None
    __state_f = None

    def __init__(self, ureg=None, atom=None, **transition):
        self._ureg = _ureg
        if ureg is not None:
            self._ureg = ureg
        if atom is not None:
            self._ureg = atom._ureg

        if "matrix_element" in transition:
            self._matrix_element = transition["matrix_element"]
        if "d" in transition:
            self._matrix_element = transition["d"]

        if "Gamma" in transition:
            Gamma = self.__ureg.Quantity(transition["Gamma"])
        elif "Aki(s^-1)" in transition:
            Gamma = self.__ureg.Quantity(float(transition["Aki(s^-1)"]), "s^-1").to(
                "Eh/ħ"
            )
        else:
            Gamma = 0.0

        if "Ei" in transition:
            Ei = self.__ureg.Quantity(transition["Ei"])

        elif "Ei(Ry)" in transition:
            Ei = self.__ureg.Quantity(
                float(remove_annotations(transition["Ei(Ry)"])), "Ry"
            ).to("Eh")
        else:
            Ei = 0.0

        if "Ef" in transition:
            Ef = self.__ureg.Quantity(transition["Ef"])
        elif "Ek(Ry)" in transition:
            Ef = self.__ureg.Quantity(
                float(remove_annotations(transition["Ek(Ry)"])), "Ry"
            ).to("Eh")
        else:
            Ef = 0.0

        super(Transition, self).__init__({"Ei": Ei, "Ef": Ef, "Gamma": Gamma})

    def __repr__(self):
        fmt = "0.4g~P"
        if self.i is not None:
            state_i = f"{self.i.valence} {self.i.term}"
        else:
            state_i = f"{self.Ei:{fmt}}"
        if self.f is not None:
            state_f = f"{self.f.valence} {self.f.term}"
        else:
            state_f = f"{self.Ef:{fmt}}"

        return (
            f"Transition({state_i} <---> {state_f}, "
            f"λ={self.λ_nm:{fmt}}, "
            f"Γ=2π×{self.Γ_MHz/(2*π):{fmt}})"
        )

    def to_dict(self):
        return {
            "i": self.__atom.states.index(self.i),
            "f": self.__atom.states.index(self.f),
            "Ei": str(self.Ei),
            "Ef": str(self.Ef),
            "Gamma": str(self.Gamma),
        }

    @property
    def _matrix_element(self):
        return self.__matrix_element

    @_matrix_element.setter
    def _matrix_element(self, matrix_element):
        if isinstance(matrix_element, str):
            matrix_element = self.__ureg(matrix_element)
        if not isinstance(matrix_element, pint.Quantity):
            matrix_element = matrix_element * self.__ureg("e a0")
        if not matrix_element.check("[charge]*[length]"):
            raise ValueError("matrix element must have units of electric dipole")
        self.__matrix_element = matrix_element

    @property
    def _Γ(self):
        ε_0 = self.__ureg.ε_0
        ħ = self.__ureg.ħ
        c = self.__ureg.c

        ω = self.ω
        J = self.f.J
        d = self.matrix_element
        return (ω ** 3) / ((3 * π * ε_0 * ħ * c ** 3) * (2 * J + 1)) * d ** 2

    @_Γ.setter
    def _Γ(self, Γ):
        ε_0 = self.__ureg.ε_0
        ħ = self.__ureg.ħ
        c = self.__ureg.c

        ω = self.ω
        J = self.f.J
        self._matrix_element = (
            (3 * π * ε_0 * ħ * c ** 3) / (ω ** 3) * (2 * J + 1) * Γ
        ) ** (1 / 2)

    @property
    def Ei(self):
        return self["Ei"]

    @property
    def Ef(self):
        return self["Ef"]

    @property
    def Gamma(self):
        return self["Gamma"]

    @property
    def Γ(self):
        return self["Gamma"]

    @property
    def Gamma_MHz(self):
        return self.Γ_MHz

    @property
    def Γ_MHz(self):
        return self.Γ.to("MHz")

    @property
    def i(self):
        return self._state_i

    @property
    def f(self):
        return self._state_f

    @property
    def ω(self):
        ħ = self.__ureg["ħ"]
        return (self.Ef - self.Ei) / ħ

    @property
    def angular_frequency(self):
        return self.ω

    @property
    def ν(self):
        return self.ω / (2 * π)

    @property
    def frequency(self):
        return self.ν

    @property
    def λ(self):
        c = self.__ureg["c"]
        try:
            return c / self.ν
        except ZeroDivisionError:
            return inf

    @property
    def wavelength(self):
        return self.λ

    @property
    def λ_nm(self):
        return self.λ.to("nm")

    @property
    def wavelength_nm(self):
        return self.λ_nm

    @property
    def saturation_intensity(self):
        h = self.__ureg["h"]
        c = self.__ureg["c"]
        return π * h * c * self.Γ / (3 * self.λ ** 3)

    @property
    def Isat(self):
        return self.saturation_intensity

    @property
    def branching_ratio(self):
        r = self.Γ * self.f.τ
        if isinstance(r, self.__ureg.Quantity):
            r = r.m
        return r

    @property
    def reduced_dipole_matrix_element(self):
        ω0 = self.ω
        ε_0 = self.__ureg["ε_0"]
        ħ = self.__ureg["ħ"]
        c = self.__ureg["c"]
        Je = self.f.J
        Γ = self.Γ
        return ((3 * π * ε_0 * ħ * c ** 3) / (ω0 ** 3) * (2 * Je + 1) * Γ) ** (1 / 2)

    @property
    def reduced_dipole_matrix_element_conjugate(self):
        Jg = self.i.J
        Je = self.f.J
        d = self.reduced_dipole_matrix_element
        return (-1) ** (Je - Jg) * d

    @property
    def d(self):
        return self.reduced_dipole_matrix_element

    @property
    def d_conj(self):
        return self.reduced_dipole_matrix_element_conjugate

    @property
    def σ0(self):
        ħ = self.__ureg["ħ"]
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
