from math import cos
from math import pi as π

import pint

from atomphys.calc.wigner import ishalfint, isint, wigner_6j


def scalar(state, omega: pint.Quantity = 0) -> pint.Quantity:
    """Calculate the scalar polarizability of a state |Ψ⟩

    Approximate the dynamical scalar polarizability α0(ω)
    of the state |Ψ⟩ by summing over matrix elements

    $$
        \\alpha^0(\\omega) =
            \\frac{1}{3(2J+1)}
            \\sum_i
                \\left(
                \\frac{1}{\\hbar (E_i - E) - \\hbar \\omega} +
                \\frac{1}{\\hbar (E_i - E) + \\hbar \\omega}
                \\right)
                \\left|\\left<\\Psi\\left|D\\right| \\Psi_i\\right>\\right|^2
    $$

    Arguments:
        state (State): state |Ψ⟩ to calculate the polarizability
        omega: Angular frequency of the field to calculate
            the dynamical polarizability. Defaults to zero,
            in which case it calculates the static polarizability.
            Must have units of `Hz`.

    Returns:
        Scalar polarizability. Has units of `C m / (V/m)`, or dipole moment / electric field.
    """
    ω = omega
    ħ = state._ureg["ħ"]
    J = state.J
    X = 1 / (3 * (2 * J + 1)) / ħ
    α = 0

    for transition in state.up:
        ω0 = transition.ω
        d = transition.reduced_dipole_matrix_element
        α += (1 / (ω0 - ω) + 1 / (ω0 + ω)) * d ** 2

    for transition in state.down:
        ω0 = -transition.ω
        d = transition.reduced_dipole_matrix_element_conjugate
        α += (1 / (ω0 - ω) + 1 / (ω0 + ω)) * d ** 2

    return (α * X).to_base_units()


def vector(state, omega: pint.Quantity = 0):
    """Calculate the vector polarizability of a state |Ψ⟩

    Approximate the dynamical vector polarizability α1(ω)
    of the state |Ψ⟩ by summing over matrix elements


    \\[
        \\alpha^1(\\omega) =
            \\sqrt{\\frac{6J}{(J+1)(2J+1)}}
            \\sum_i
                \\left\\{
                \\begin{matrix}
                1 & 1 & 1 \\\\ J & J & J_i
                \\end{matrix}
                \\right\\}
                (-1)^{1+J+J_i}
                \\left(
                \\frac{1}{\\hbar (E_i - E) - \\hbar \\omega} +
                \\frac{1}{\\hbar (E_i - E) + \\hbar \\omega}
                \\right)
                \\left|\\left<\\Psi\\left|D\\right| \\Psi_i\\right>\\right|^2
    \\]

    Arguments:
        state (State): state |Ψ⟩ to calculate the polarizability
        omega: Angular frequency of the field to calculate
            the dynamical polarizability. Defaults to zero,
            in which case it calculates the static polarizability.
            Must have units of `Hz`.

    Returns:
        Vector polarizability. Has units of `C m / (V/m)`, or dipole moment / electric field.
    """
    ω = omega
    ħ = state._ureg["ħ"]
    J = state.J
    X = ((6 * J) / (4 * (2 * J + 1) * (J + 1))) ** (1 / 2) / ħ
    α = 0

    for transition in state.up:
        ω0 = transition.ω
        d = transition.reduced_dipole_matrix_element
        Jp = transition.f.J
        sixJ = wigner_6j(1, 1, 1, J, J, Jp)
        α += (-1) ** (J + Jp + 1) * sixJ * (1 / (ω0 - ω) + 1 / (ω0 + ω)) * d ** 2

    for transition in state.down:
        ω0 = -transition.ω
        d = transition.reduced_dipole_matrix_element_conjugate
        Jp = transition.i.J
        sixJ = wigner_6j(1, 1, 1, J, J, Jp)
        α += -((-1) ** (J + Jp + 1)) * sixJ * (1 / (ω0 - ω) + 1 / (ω0 + ω)) * d ** 2

    return (α * X).to_base_units()


def tensor(state, omega: pint.Quantity = 0):
    """Calculate the tensor polarizability of a state |Ψ⟩

    Approximate the dynamical tensor polarizability α0(ω)
    of the state |Ψ⟩ by summing over matrix elements

    \\[
        \\alpha^2(\\omega) = -\\sqrt{\\frac{20J(2J-1)}{6(J+1)(2J+1)(2J+3)}}
            \\sum_i
                \\left\\{
                \\begin{matrix}
                1 & 1 & 2 \\\\ J & J & J_i
                \\end{matrix}
                \\right\\}
                (-1)^{1+J+J_i}
                \\left(
                \\frac{1}{\\hbar (E_i - E) - \\hbar \\omega} +
                \\frac{1}{\\hbar (E_i - E) + \\hbar \\omega}\\right)
                \\left|
                \\left<\\Psi\\left|D\\right| \\Psi_i\\right>\\right|^2
    \\]

    Arguments:
        state (State): state |Ψ⟩ to calculate the polarizability
        omega: Angular frequency of the field to calculate
            the dynamical polarizability. Defaults to zero,
            in which case it calculates the static polarizability.
            Must have units of `Hz`.

    Returns:
        Tensor polarizability. Has units of `C m / (V/m)`, or dipole moment / electric field.
    """
    ω = omega
    ħ = state._ureg["ħ"]
    J = state.J
    X = (
        -(
            ((20 * J * (2 * J - 1)) / (6 * (J + 1) * (2 * J + 1) * (2 * J + 3)))
            ** (1 / 2)
        )
        / ħ
    )
    α = 0

    for transition in state.up:
        ω0 = transition.ω
        d = transition.reduced_dipole_matrix_element
        Jp = transition.f.J
        sixJ = wigner_6j(1, 1, 2, J, J, Jp)
        α += (-1) ** (J + Jp + 1) * sixJ * (1 / (ω0 - ω) + 1 / (ω0 + ω)) * d ** 2

    for transition in state.down:
        ω0 = -transition.ω
        d = transition.reduced_dipole_matrix_element_conjugate
        Jp = transition.i.J
        sixJ = wigner_6j(1, 1, 2, J, J, Jp)
        α += (-1) ** (J + Jp + 1) * sixJ * (1 / (ω0 - ω) + 1 / (ω0 + ω)) * d ** 2

    return (α * X).to_base_units()


def total(
    state,
    mJ: float = None,
    omega: pint.Quantity = 0,
    A: float = 0,
    theta_k: float = 0,
    theta_p: float = π / 2,
) -> pint.Quantity:
    """Calculate the polarizability of a state for a given field polarization

    Calculates the dynamical polarizability α(ω) by taking the sum of the
    scalar, vector, and tensor poarts according to

    \\[
        \\alpha(\\omega) = \\alpha^0(\\omega) +
        A \\cos(\\theta_k) \\frac{m_J}{J} \\alpha^1(\\omega) +
        \\frac{1}{2}\\left(3 \\cos^2(\\theta_p) - 1\\right)\\frac{3m_J^2 - J(J+1)}{J(2J-1)}\\alpha^2(\\omega)
    \\]

    Arguments:
        state (State): state |Ψ⟩ to calculate the polarizability
        mJ: Zeeman sublevel to calculate the vector and tensor
            polarizability. If `None` only calculates the scalar
            component. Must be an integer if J is an integer, or
            a half-integer if J is a half integer.
        omega: angular frequency of the field to calculate
            the dynamical polarizability. Defaults to zero,
            in which case it calculates the static polarizability.
            Must have units of `Hz`.
        A: degree of circular polarization, with ±1 being circular
            polarization and 0 being linear polarization
        theta_k: angle between wave vector and quantization axis
        theta_p: angle between polarization vector and quantization axis

    Returns:
        Polarizability. Has units of `C m / (V/m)`, or dipole moment / electric field.

    Raises:
        ZeroDivisionError
    """
    θk = theta_k
    θp = theta_p

    J = state.J
    α0 = scalar(state, omega)
    α1 = vector(state, omega)
    α2 = tensor(state, omega)
    if mJ is None:
        return α0
    if not (isint(mJ) or (ishalfint(mJ) and ishalfint(J))):
        raise ValueError(
            "mJ must be an integer, or a half-integer if J is a half-integer"
        )

    C1 = A * cos(θk) * mJ / J

    try:
        C2 = (
            (3 * cos(θp) ** 2 - 1) / 2 * (3 * mJ ** 2 - J * (J + 1)) / (J * (2 * J - 1))
        )
    except ZeroDivisionError:
        C2 = 0

    return (α0 + C1 * α1 + C2 * α2).to_base_units()
