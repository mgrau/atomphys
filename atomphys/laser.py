from math import inf
from math import pi as π

import pint

from . import _ureg
from .util import default_units


class Laser:
    _ureg: pint.UnitRegistry
    __omega = pint.Quantity
    __linewidth = pint.Quantity
    __electric_field = pint.Quantity

    __A = 0
    __theta_k = 0
    __theta_p = π / 2

    def __init__(self, ureg=None, laser=None, **new_laser):
        if ureg is not None:
            self._ureg = ureg
        else:
            self._ureg = _ureg

        self.omega = 0
        self.__linewidth = 0
        self.__electric_field = 0

        if laser is not None:
            self._ureg = laser._ureg
            self.__omega = laser.__omega
            self.__linewidth = laser.__linewidth
            self.__electric_field = laser.__electric_field

            self.__A = laser.__A
            self.__theta_k = laser.__theta_k
            self.__theta_p = laser.__theta_p

        for attr in new_laser:
            if attr in dir(self):
                self.__setattr__(attr, new_laser[attr])

    def __repr__(self):
        fmt = "0.4g~P"
        return f"Laser(λ={self.λ:{fmt}})"

    # ---------
    # Frequency
    # ---------

    @property
    def omega(self):
        return self.__omega

    @omega.setter
    @default_units("THz")
    def omega(self, ω):
        self.__omega = ω

    @property
    def angular_frequency(self):
        return self.omega

    @angular_frequency.setter
    def angular_frequency(self, ω):
        self.omega = ω

    @property
    def ω(self):
        return self.omega

    @ω.setter
    def ω(self, ω):
        self.omega = ω

    @property
    def ν(self):
        return self.omega / (2 * π)

    @ν.setter
    def ν(self, ν):
        self.omega = 2 * π * ν

    @property
    def nu(self):
        return self.ν

    @nu.setter
    def nu(self, ν):
        self.ν = ν

    @property
    def frequency(self):
        return self.ν

    @frequency.setter
    def frequency(self, ν):
        self.ν = ν

    # ----------
    # Wavelength
    # ----------

    @property
    def wavelength(self):
        c = self._ureg.c
        try:
            return (2 * π * c / self.omega).to("nm")
        except ZeroDivisionError:
            return inf * self._ureg("nm")

    @wavelength.setter
    @default_units("nm")
    def wavelength(self, λ):
        c = self._ureg.c
        self.omega = 2 * π * c / λ

    @property
    def λ(self):
        return self.wavelength

    @λ.setter
    def λ(self, λ):
        self.wavelength = λ

    # ---------
    # Linewidth
    # ---------

    @property
    def linewidth(self):
        return self.__linewidth

    @linewidth.setter
    @default_units("Hz")
    def linewidth(self, linewidth):
        self.__linewidth = linewidth

    # --------------
    # Electric Field
    # --------------

    @property
    def electric_field(self):
        return self.__electric_field

    @electric_field.setter
    @default_units("V/m")
    def electric_field(self, E):
        self.__electric_field = E

    @property
    def E(self):
        return self.electric_field

    @E.setter
    def E(self, E):
        self.electric_field = E

    @property
    def intensity(self):
        c = self._ureg.c
        ε_0 = self._ureg.ε_0
        return self.electric_field ** 2 * (c * ε_0 / 2)

    @intensity.setter
    @default_units("W/m^2")
    def intensity(self, I):
        c = self._ureg.c
        ε_0 = self._ureg.ε_0
        self.electric_field = (2 * I / (c * ε_0)) ** (1 / 2)

    @property
    def I(self):
        return self.intensity

    @I.setter
    def I(self, I):
        self.intensity = I

    # ------------
    # Polarization
    # ------------

    @property
    def A(self):
        return self.__A

    @A.setter
    def A(self, A):
        self.__A = A

    @property
    def theta_k(self):
        return self.__theta_k

    @theta_k.setter
    def theta_k(self, theta_k):
        self.__theta_k = theta_k

    @property
    def theta_p(self):
        return self.__theta_p

    @theta_p.setter
    def theta_p(self, theta_p):
        self.__theta_p = theta_p

    # ---------------------
    # high level properties
    # ---------------------

    def Rabi_frequency(self, transition, Rabi_frequency=None):
        # this is not quite right, as d is reduced dipole matrix element.
        # also transition is not necessarily dipole
        ħ = self._ureg.ħ
        if Rabi_frequency is None:
            return (transition.d * self.E / ħ).to_base_units()
        else:
            if not Rabi_frequency.check("1/[time]"):
                raise ValueError("Rabi_frequency must be a frequency")
            self.E = ħ * Rabi_frequency / transition.d

    @property
    def Ω(self):
        return self.Rabi_frequency
