
from . import _ureg as _units
from math import pi as π
from math import inf


class Laser:
    __units = {}
    __omega = None
    __linewidth = None

    __A = 0
    __theta_k = 0
    __theta_p = π/2
    __electric_field = None

    def __init__(self, units=None, laser=None, **new_laser):
        if units:
            self.__units = units
        else:
            self.__units = _units

        self.__omega = self.__units('0 Hz')
        self.__linewidth = self.__units('0 Hz')
        self.__intensity = self.__units('0 V/m')

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
        fmt = '0.4g~P'
        return (
            f'Laser(λ={self.λ:{fmt}})'
        )

    @property
    def ω(self):
        return self.__omega

    @ω.setter
    def ω(self, ω):
        if not ω.check('[frequency]'):
            raise ValueError('ω must be a frequency')
        self.__omega = ω

    @property
    def omega(self):
        return self.ω

    @omega.setter
    def omega(self, ω):
        if not ω.check('[frequency]'):
            raise ValueError('omega must be a frequency')
        self.ω = ω

    @property
    def angular_frequency(self):
        return self.ω

    @angular_frequency.setter
    def angular_frequency(self, ω):
        if not ω.check('[frequency]'):
            raise ValueError('angular frequency must be a frequency')
        self.ω = ω

    @property
    def ν(self):
        return self.ω/(2*π)

    @ν.setter
    def ν(self, ν):
        if not ν.check('[frequency]'):
            raise ValueError('ν must be a frequency')
        self.ω = 2*π*ν

    @property
    def nu(self):
        return self.ν

    @nu.setter
    def nu(self, ν):
        if not ν.check('[frequency]'):
            raise ValueError('nu must be a frequency')
        self.ν = ν

    @property
    def frequency(self):
        return self.ν

    @frequency.setter
    def frequency(self, ν):
        if not ν.check('[frequency]'):
            raise ValueError('frequency must be a frequency')
        self.ν = ν

    @property
    def λ(self):
        c = self.__units['c']
        try:
            return c/self.ν
        except ZeroDivisionError:
            return inf*self.ν

    @λ.setter
    def λ(self, λ):
        if not λ.check('[length]'):
            raise ValueError('λ must be a length')
        c = self.__units['c']
        self.ν = c/λ

    @property
    def wavelength(self):
        return self.λ

    @wavelength.setter
    def wavelength(self, λ):
        if not λ.check('[length]'):
            raise ValueError('wavelength must be a length')
        self.λ = λ

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

    @property
    def E(self):
        return self.__electric_field

    @E.setter
    def E(self, E):
        if not E.check('[mass]*[length]/[current]/[time]^3'):
            raise ValueError('E must be an electric field')
        self.__electric_field = E

    @property
    def electric_field(self):
        return self.E

    @electric_field.setter
    def electric_field(self, E):
        if not E.check('[mass]*[length]/[current]/[time]^3'):
            raise ValueError('electric_field must be an electric field')
        self.E = E

    @property
    def I(self):
        c = self.__units['c']
        ε_0 = self.__units['ε_0']
        return self.E**2*(c*ε_0/2)

    @I.setter
    def I(self, I):
        if not I.check('[mass]/[time]^3'):
            raise ValueError('I must be an intensity')
        c = self.__units['c']
        ε_0 = self.__units['ε_0']
        self.E = (2*I/(c*ε_0))**(1/2)

    @property
    def intensity(self):
        return self.I

    @intensity.setter
    def intensity(self, I):
        if not I.check('[mass]/[time]^3'):
            raise ValueError('intensity must be an intensity')
        self.I = I

    def Ω(self, transition, Rabi_frequency=None):
        # this is not quite right, as d is reduced dipole matrix element.
        # also transition is not necessarily dipole
        ħ = self.__units['ħ']
        if Rabi_frequency is None:
            return (transition.d * self.E / ħ).to_base_units()
        else:
            if not Rabi_frequency.check('1/[time]'):
                raise ValueError('Rabi_frequency must be a frequency')
            self.E = ħ*Rabi_frequency/transition.d

    @property
    def Rabi_frequency(self):
        return self.Ω
