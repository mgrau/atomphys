from math import pi as π

from .wigner import wigner_6j


def scalar(state, omega=0):
    ω = omega
    ħ = state._ureg['ħ']
    α = 0
    for transition in state.up:
        ω0 = transition.ω
        d = transition.reduced_dipole_matrix_element
        α += (2/3)*ω0/(ω0**2-ω**2)*d**2/ħ

    for transition in state.down:
        ω0 = -transition.ω
        d = transition.reduced_dipole_matrix_element_conjugate
        α += (2/3)*ω0/(ω0**2-ω**2)*d**2/ħ

    return α


def vector(state, omega=0):
    ω = omega
    ħ = state._ureg['hbar']
    α = 0
    J = state.J
    X = ((6*J*(2*J+1)/(J+1)))**(1/2)

    for transition in state.up:
        ω0 = transition.ω
        d = transition.reduced_dipole_matrix_element
        Jp = transition.f.J
        sixJ = wigner_6j(1, 1, 1, J, J, Jp)
        α += (-1)**(J+Jp+1)*X*sixJ*ω0/(ω0**2-ω**2)*d**2/ħ

    for transition in state.down:
        ω0 = -transition.ω
        d = transition.reduced_dipole_matrix_element_conjugate
        Jp = transition.i.J
        sixJ = wigner_6j(1, 1, 1, J, J, Jp)
        α += (-1)**(J+Jp+1)*X*sixJ*ω0/(ω0**2-ω**2)*d**2/ħ

    return α


def tensor(state, omega=0):
    ω = omega
    ħ = state._ureg['hbar']
    α = 0
    J = state.J
    X = ((40*J*(2*J+1)*(2*J-1))/(3*(J+1)*(2*J+3)))**(1/2)

    for transition in state.up:
        ω0 = transition.ω
        d = transition.reduced_dipole_matrix_element
        Jp = transition.f.J
        sixJ = wigner_6j(1, 1, 2, J, J, Jp)
        α += (-1)**(J+Jp+1)*X*sixJ*ω0/(ω0**2-ω**2)*d**2/ħ

    for transition in state.down:
        ω0 = -transition.ω
        d = transition.reduced_dipole_matrix_element_conjugate
        Jp = transition.i.J
        sixJ = wigner_6j(1, 1, 2, J, J, Jp)
        α += (-1)**(J+Jp+1)*X*sixJ*ω0/(ω0**2-ω**2)*d**2/ħ

    return α
