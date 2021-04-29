
from . import _ureg as _units
from math import pi as π


class Laser:
    __units = {}
    __omega = None

    def __init__(self, units=None, **laser):
        if units:
            self.__units = units
        else:
            self.__units = _units

        for attr in laser:
            if attr in dir(self):
                self.__setattr__(attr, laser[attr])

    @property
    def ω(self):
        return self.__omega

    @ω.setter
    def ω(self, ω):
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
        return self.ω/(2*π)

    @ν.setter
    def ν(self, ν):
        self.ω = 2*π*ν

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

    @property
    def λ(self):
        c = self.__units['c']
        try:
            return c/self.ν
        except ZeroDivisionError:
            return inf*self.ν

    @λ.setter
    def λ(self, λ):
        c = self.__units['c']
        self.ν = c/λ

    @property
    def wavelength(self):
        return self.λ

    @wavelength.setter
    def wavelength(self, λ):
        self.λ = λ
