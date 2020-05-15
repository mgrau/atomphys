from .data import State, Transition, get_states, get_transitions
from math import pi as π
import numpy as np


try:
    from . import _ureg, _Q, _HAS_PINT
except ImportError:
    _HAS_PINT = False
    _ureg = None
    _Q = None


class Atom():

    def __init__(self, atom, USE_UNITS=True):
        self.USE_UNITS = USE_UNITS
        self._states = get_states(atom, USE_UNITS=USE_UNITS)
        self._transitions = get_transitions(atom, USE_UNITS=USE_UNITS)

        # sort by upper state energy first
        self._transitions.sort(key=lambda transition: transition.Ef)
        # sort then by lower state energy
        # because sort is stable, this produces a list sorted by both upper and lower state energy
        self._transitions.sort(key=lambda transition: transition.Ei)

        # index the transitions to the states
        for transition in self._transitions:
            transition['i'] = next(
                state for state in self._states if state.energy == transition.Ei)
            transition['f'] = next(
                state for state in self._states if state.energy == transition.Ef)

    @property
    def states(self):
        return self._states

    @property
    def transitions(self):
        return self._transitions

    def transitions_down(self, state):
        if not isinstance(state, State):
            # treat state like an index instead
            state = self._states[state]

        transitions = [
            transition for transition in self._transitions if transition.Ef == state.energy]
        transitions.sort(key=lambda transition: transition.Ei)
        return transitions

    def transitions_up(self, state):
        if not isinstance(state, State):
            # treat state like an index instead
            state = self._states[state]

        transitions = [
            transition for transition in self._transitions if transition.Ei == state.energy]
        transitions.sort(key=lambda transition: transition.Ef)
        return transitions

    def lifetime(self, state):
        transitions = self.transitions_down(state)
        if not transitions:
            return -1
        else:
            Gamma = [transition.Gamma for transition in transitions]
            return 1/sum(Gamma)

    def scalar_polarizability(self, state, ω=0):
        if not isinstance(state, State):
            # treat state like an index instead
            state = self._states[state]

        if not hasattr(ω, '__len__'):
            ω = [ω]

        if _HAS_PINT and self.USE_UNITS:
            hbar = _ureg['hbar']
            ε_0 = _ureg['ε_0']
            c = _ureg['c']
        else:
            hbar = 1
            ε_0 = 1
            c = 1

        transitions_up = self.transitions_up(state)
        ω0 = [(t.Ef - t.Ei)/hbar for t in transitions_up]
        Γ = [t.Gamma for t in transitions_up]
        deg = np.array([(2*t.f.J + 1)/(2*t.i.J + 1)
                        for t in transitions_up])
        if _HAS_PINT and self.USE_UNITS:
            if len(ω0):
                ω0 = _Q(np.array([ω.magnitude for ω in ω0]), ω0[0].units)
            else:
                ω0 = _Q([], _ureg['Eh/hbar'])

            if len(Γ):
                Γ = _Q(np.array([γ.magnitude for γ in Γ]), Γ[0].units)
            else:
                Γ = _Q([], _ureg['Eh/hbar'])
        else:
            ω0 = np.array(ω0)
            Γ = np.array(Γ)

        RME2 = 3*π*ε_0*c**3 * ω0**-3 * deg * Γ

        RME2 = np.tile(RME2, (np.size(ω), 1))
        ω0_m = np.tile(ω0, (np.size(ω), 1))
        ω_m = np.tile(ω, (np.size(ω0), 1)).T
        Δ = ω0_m/(ω0_m**2 - ω_m**2)
        α0_up = np.sum((2/3)*RME2*Δ, axis=1)

        transitions_down = self.transitions_down(state)
        ω0 = [(t.Ef - t.Ei)/hbar for t in transitions_down]
        Γ = [t.Gamma for t in transitions_down]
        if _HAS_PINT and self.USE_UNITS:
            if len(ω0):
                ω0 = _Q(np.array([ω.magnitude for ω in ω0]), ω0[0].units)
            else:
                ω0 = _Q([], _ureg['Eh/hbar'])

            if len(Γ):
                Γ = _Q(np.array([γ.magnitude for γ in Γ]), Γ[0].units)
            else:
                Γ = _Q([], _ureg['Eh/hbar'])
        else:
            ω0 = np.array(ω0)
            Γ = np.array(Γ)

        RME2 = 3*π*ε_0*c**3 * ω0**-3 * Γ

        RME2 = np.tile(RME2, (np.size(ω), 1))
        ω0_m = np.tile(ω0, (np.size(ω), 1))
        ω_m = np.tile(ω, (np.size(ω0), 1)).T
        Δ = ω0_m/(ω0_m**2 - ω_m**2)
        α0_down = np.sum((2/3)*RME2*Δ, axis=1)

        α0 = α0_up - α0_down
        return α0
