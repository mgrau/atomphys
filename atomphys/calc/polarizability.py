from math import pi as π


def scalar(state, omega=0):
    ω = omega

    ε_0 = state._ureg['ε_0']
    c = state._ureg['c']

    α0 = 0
    for transition in state.transitions:
        ω0 = transition.ω
        Γ = transition.Γ
        if transition in state.up:
            deg = (2*transition.f.J + 1)/(2*transition.i.J + 1)
        else:
            deg = -1
        RME2 = 3*π*ε_0*c**3 * ω0**-3 * deg * Γ
        Δ = ω0/(ω0**2 - ω**2)
        α0 += (2/3)*RME2*Δ

    return α0
