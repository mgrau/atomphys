from math import inf
from math import pi as π

from . import _ureg as _units


class Laser:
    __units = {}
    __omega = None
    __linewidth = None

    __A = 0
    __theta_k = 0
    __theta_p = π / 2
    __electric_field = None

    def __init__(self, units=None, laser=None, **new_laser):
        if units:
            self.__units = units
        else:
            self.__units = _units

        self.__omega = self.__units("0 Hz")
        self.__linewidth = self.__units("0 Hz")
        self.__intensity = self.__units("0 V/m")

        if laser is not None:
            self.__units = laser.__units
            self.__omega = laser.__omega
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
    def ω(self):
        return self.__omega

    @ω.setter
    def ω(self, ω):
        if isinstance(ω, str):
            ω = self.__units(ω)
        if not isinstance(ω, self.__units.Quantity):
            ω = ω * self.__units.Hz
        if not ω.check("[frequency]"):
            raise ValueError("ω must have units of frequency")
        self.__omega = ω

    @property
    def omega(self):
        return self.ω

    @omega.setter
    def omega(self, ω):
        self.ω = ω

    @property
    def angular_frequency(self):
        return self.ω

    @angular_frequency.setter
    def angular_frequency(self, ω):
        self.ω = ω

    @property
    def ν(self):
        return self.ω / (2 * π)

    @ν.setter
    def ν(self, ν):
        self.ω = 2 * π * ν

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
    def λ(self):
        c = self.__units["c"]
        try:
            return c / self.ν
        except ZeroDivisionError:
            return inf * self.ν

    @λ.setter
    def λ(self, λ):
        if isinstance(λ, str):
            λ = self.__units(λ)
        if not isinstance(λ, self.__units.Quantity):
            λ = λ * self.__units.nm
        if not λ.check("[length]"):
            raise ValueError("wavelength must have units of length")
        c = self.__units["c"]
        self.ν = c / λ

    @property
    def wavelength(self):
        return self.λ

    @wavelength.setter
    def wavelength(self, λ):
        self.λ = λ

    # --------------
    # Electric Field
    # --------------

    @property
    def E(self):
        return self.__electric_field

    @E.setter
    def E(self, E):
        if isinstance(E, str):
            E = self.__units(E)
        if not isinstance(E, self.__units.Quantity):
            E = E * self.__units("V/m")
        if not E.check("[mass]*[length]/[current]/[time]^3"):
            raise ValueError("E must have units of electric field")
        self.__electric_field = E

    @property
    def electric_field(self):
        return self.E

    @electric_field.setter
    def electric_field(self, E):
        self.E = E

    @property
    def I(self):
        c = self.__units["c"]
        ε_0 = self.__units["ε_0"]
        return self.E ** 2 * (c * ε_0 / 2)

    @I.setter
    def I(self, I):
        if isinstance(I, str):
            I = self._units(I)
        if not isinstance(I, self.__units.Quantity):
            I = I * self.__units("W/m^2")
        if not I.check("[mass]/[time]^3"):
            raise ValueError("I must have units of intensity")
        c = self.__units["c"]
        ε_0 = self.__units["ε_0"]
        self.E = (2 * I / (c * ε_0)) ** (1 / 2)

    @property
    def intensity(self):
        return self.I

    @intensity.setter
    def intensity(self, I):
        self.I = I

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

    def Ω(self, transition, Rabi_frequency=None):
        # this is not quite right, as d is reduced dipole matrix element.
        # also transition is not necessarily dipole
        ħ = self.__units["ħ"]
        if Rabi_frequency is None:
            return (transition.d * self.E / ħ).to_base_units()
        else:
            if not Rabi_frequency.check("1/[time]"):
                raise ValueError("Rabi_frequency must be a frequency")
            self.E = ħ * Rabi_frequency / transition.d

    @property
    def Rabi_frequency(self):
        return self.Ω
