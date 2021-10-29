from math import inf

import pint
import pytest

from atomphys import _ureg
from atomphys.laser import Laser
from atomphys.state import State
from atomphys.transition import Transition


def test_initialization():
    assert Laser(omega=100).omega == _ureg("100 THz")
    assert Laser(omega=200).ω == _ureg("200 THz")
    assert Laser(omega=300).angular_frequency == _ureg("300 THz")
    assert Laser(ω="100 GHz").omega == _ureg("100 GHz")
    assert Laser(angular_frequency=_ureg("1 Hz")).omega == _ureg("1 Hz")
    assert Laser(ν=100).λ == _ureg.Quantity(2997.92458, "nm")
    assert Laser(nu=1000).wavelength == _ureg.Quantity(299.792458, "nm")
    assert Laser(frequency=10).wavelength == _ureg.Quantity(29979.2458, "nm")
    assert Laser(nu=1).ν == _ureg("1 THz")
    assert Laser(ν=2).nu == _ureg("2 THz")
    assert Laser(frequency=3).frequency == _ureg("3 THz")

    assert Laser(ω=0).wavelength == _ureg.Quantity(inf, "nm")

    assert Laser(linewidth=1000).linewidth == _ureg.Quantity(1, "kHz")

    assert Laser(electric_field=1).electric_field == _ureg.Quantity(1, "V/m")
    assert Laser(electric_field="2 V/cm").E == _ureg.Quantity(2, "V/cm")
    assert Laser(E="3 V/m").electric_field == _ureg.Quantity(3, "V/m")
    assert Laser(E=1).intensity == _ureg.Quantity(0.5, "c*ε_0*(V/m)^2")
    assert Laser(electric_field="2 V/m").I == _ureg.Quantity(2, "c*ε_0*(V/m)^2")

    assert Laser(A=1).A == 1
    assert Laser(theta_k=2).theta_k == 2
    assert Laser(theta_p=3).theta_p == 3

    assert Laser(laser=Laser(ω=400), omega=500).angular_frequency == _ureg("500 THz")


def test_repr():
    assert str(Laser(λ=532)) == "Laser(λ=532 nm)"
    assert str(Laser(nu=1000)) == "Laser(λ=299.8 nm)"


def test_units():
    ureg = pint.UnitRegistry()
    assert Laser(λ=500).λ - Laser(λ=400).λ == _ureg("100 nm")

    with pytest.raises(ValueError):
        Laser(λ=500).λ - Laser(λ=400, ureg=ureg).λ

    assert Laser()._ureg is _ureg
    assert Laser()._ureg is not ureg
    assert Laser(ureg=ureg)._ureg is not _ureg
    assert Laser(ureg=ureg)._ureg is ureg


def test_Rabi_frequency():
    transition = Transition(
        state_i=State("2S1/2"),
        state_f=State("2P1/2", energy="h*c/(532 nm)"),
        d="1 e a0",
    )
    laser = Laser(I="1 mW/cm^2")
    laser.Rabi_frequency(transition).m_as("MHz") == pytest.approx(6.978557)
    laser.Ω(transition).m_as("MHz") == pytest.approx(6.978557)

    laser.set_Rabi_frequency("100 MHz", transition)
    assert laser.I.m_as("mW/cm^2") == pytest.approx(205.337706)
    laser.set_Ω("10 MHz", transition)
    assert laser.intensity.m_as("mW/cm^2") == pytest.approx(2.05337706)
