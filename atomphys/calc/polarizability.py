from math import pi as π
from math import cos

from .wigner import wigner_6j


def scalar(state, omega=0):
    ω = omega
    ħ = state._ureg['ħ']
    J = state.J
    X = 1/(3*(2*J+1))/ħ
    α = 0

    for transition in state.up:
        ω0 = transition.ω
        d = transition.reduced_dipole_matrix_element
        α += (1/(ω0-ω) + 1/(ω0+ω))*d**2

    for transition in state.down:
        ω0 = -transition.ω
        d = transition.reduced_dipole_matrix_element_conjugate
        α += (1/(ω0-ω) + 1/(ω0+ω))*d**2

    return (α*X).to_base_units()


def vector(state, omega=0):
    ω = omega
    ħ = state._ureg['ħ']
    J = state.J
    X = ((6*J)/(4*(2*J+1)*(J+1)))**(1/2)/ħ
    α = 0

    for transition in state.up:
        ω0 = transition.ω
        d = transition.reduced_dipole_matrix_element
        Jp = transition.f.J
        sixJ = wigner_6j(1, 1, 1, J, J, Jp)
        α += (-1)**(J+Jp+1)*sixJ*(1/(ω0-ω) + 1/(ω0+ω))*d**2

    for transition in state.down:
        ω0 = -transition.ω
        d = transition.reduced_dipole_matrix_element_conjugate
        Jp = transition.i.J
        sixJ = wigner_6j(1, 1, 1, J, J, Jp)
        α += -(-1)**(J+Jp+1)*sixJ*(1/(ω0-ω) + 1/(ω0+ω))*d**2

    return (α*X).to_base_units()


def tensor(state, omega=0):
    ω = omega
    ħ = state._ureg['ħ']
    J = state.J
    X = -((20*J*(2*J-1))/(6*(J+1)*(2*J+1)*(2*J+3)))**(1/2)/ħ
    α = 0

    for transition in state.up:
        ω0 = transition.ω
        d = transition.reduced_dipole_matrix_element
        Jp = transition.f.J
        sixJ = wigner_6j(1, 1, 2, J, J, Jp)
        α += (-1)**(J+Jp+1)*sixJ*(1/(ω0-ω) + 1/(ω0+ω))*d**2

    for transition in state.down:
        ω0 = -transition.ω
        d = transition.reduced_dipole_matrix_element_conjugate
        Jp = transition.i.J
        sixJ = wigner_6j(1, 1, 2, J, J, Jp)
        α += (-1)**(J+Jp+1)*sixJ*(1/(ω0-ω) + 1/(ω0+ω))*d**2

    return (α*X).to_base_units()


def total(state, mJ=None, omega=0, A=0, theta_k=0, theta_p=π/2):
    θk = theta_k
    θp = theta_p

    J = state.J
    α0 = scalar(state, omega)
    α1 = vector(state, omega)
    α2 = tensor(state, omega)
    if mJ is None:
        return α0

    C1 = A*cos(θk)*mJ/J

    try:
        C2 = (3*cos(θp)**2-1)/2 * (3*mJ**2 - J*(J+1)) / (J*(2*J-1))
    except ZeroDivisionError:
        C2 = 0

    return (α0 + C1*α1 + C2*α2).to_base_units()
