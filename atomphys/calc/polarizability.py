from math import pi as π

from sympy import N
from sympy.physics.wigner import wigner_6j


def scalar(state, omega=0):
    ω = omega
    ħ = state._ureg['ħ']
    α0 = 0
    for transition in state.up:
        ω0 = transition.ω
        d = transition.reduced_dipole_matrix_element
        α0 += (2/3)*ω0/(ω0**2-ω**2)*d**2/ħ

    for transition in state.down:
        ω0 = -transition.ω
        d = transition.reduced_dipole_matrix_element_conjugate
        α0 += (2/3)*ω0/(ω0**2-ω**2)*d**2/ħ

    return α0


def vector(state, omega=0):
    ω = omega
    ħ = state._ureg['hbar']
    α0 = 0
    J = state.J
    X = ((6*J*(2*J+1)/(J+1)))**(1/2)

    for transition in state.up:
        ω0 = transition.ω
        d = transition.reduced_dipole_matrix_element
        Jp = transition.f.J
        sixJ = N(wigner_6j(1, 1, 1, J, J, Jp))
        α0 += (-1)**(J+Jp+1)*X*sixJ*ω0/(ω0**2-ω**2)*d**2/ħ

    for transition in state.down:
        ω0 = -transition.ω
        d = transition.reduced_dipole_matrix_element_conjugate
        Jp = transition.i.J
        sixJ = N(wigner_6j(1, 1, 1, J, J, Jp))
        α0 += (-1)**(J+Jp+1)*X*sixJ*ω0/(ω0**2-ω**2)*d**2/ħ

    return α0


def tensor(state, omega=0):
    ω = omega
    ħ = state._ureg['hbar']
    α0 = 0
    J = state.J
    X = ((40*J*(2*J+1)*(2*J-1))/(3*(J+1)*(2*J+3)))**(1/2)

    for transition in state.up:
        ω0 = transition.ω
        d = transition.reduced_dipole_matrix_element
        Jp = transition.f.J
        sixJ = N(wigner_6j(1, 1, 2, J, J, Jp))
        α0 += (-1)**(J+Jp+1)*X*sixJ*ω0/(ω0**2-ω**2)*d**2/ħ

    for transition in state.down:
        ω0 = -transition.ω
        d = transition.reduced_dipole_matrix_element_conjugate
        Jp = transition.i.J
        sixJ = N(wigner_6j(1, 1, 2, J, J, Jp))
        α0 += (-1)**(J+Jp+1)*X*sixJ*ω0/(ω0**2-ω**2)*d**2/ħ

    return α0
