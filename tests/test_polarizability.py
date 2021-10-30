import pytest

from atomphys import State, Transition
from atomphys.calc import polarizability


def test_scalar():
    state = State("2P1/2", energy="2eV")
    Transition(State("2S1/2", energy="0 eV"), state, d="1 e a0")
    Transition(state, State("2S1/2", energy="4 eV"), d="1 e a0")
    assert polarizability.scalar(state) == 0
    assert polarizability.scalar(state, omega=state._ureg.Quantity(2000, "THz")) == 0
    assert polarizability.scalar(state, omega=state._ureg.Quantity(4000, "THz")) == 0

    state = State("2P1/2", energy="2eV")
    Transition(State("2S1/2", energy="0 eV"), state, d="1 e a0")
    assert polarizability.scalar(state).m_as("(e a0)^2/E_h") == pytest.approx(
        -4.5352310
    )
    assert polarizability.scalar(state, omega=state._ureg.Quantity(2000, "THz")).m_as(
        "(e a0)^2/E_h"
    ) == pytest.approx(-8.002073)
    assert polarizability.scalar(state, omega=state._ureg.Quantity(4000, "THz")).m_as(
        "(e a0)^2/E_h"
    ) == pytest.approx(6.187455)


def test_vector():
    state = State("3P0", energy="2eV")
    Transition(State("1S0", energy="0 eV"), state, d="1 e a0")
    assert polarizability.vector(state) == 0

    state = State("3P1", energy="2eV")
    Transition(State("1S0", energy="0 eV"), state, d="1 e a0")
    assert polarizability.vector(state).m_as("(e a0)^2/E_h") == pytest.approx(-4.535231)

    state = State("2P1/2", energy="2eV")
    Transition(State("2S1/2", energy="0 eV"), state, d="1 e a0")
    Transition(state, State("2S1/2", energy="4 eV"), d="1 e a0")
    assert polarizability.vector(state).m_as("(e a0)^2/E_h") == pytest.approx(-9.070462)


def test_tensor():
    state = State("2P1/2", energy="2eV")
    Transition(State("2S1/2", energy="0 eV"), state, d="1 e a0")
    Transition(state, State("2S1/2", energy="5 eV"), d="1 e a0")
    assert polarizability.tensor(state) == 0

    state = State("2P3/2", energy="2eV")
    Transition(State("2S1/2", energy="0 eV"), state, d="1 e a0")
    Transition(state, State("2S1/2", energy="4 eV"), d="1 e a0")
    assert polarizability.tensor(state) == 0

    state = State("2P3/2", energy="2eV")
    Transition(State("2S1/2", energy="0 eV"), state, d="1 e a0")
    Transition(state, State("2S1/2", energy="5 eV"), d="1 e a0")
    assert polarizability.tensor(state).m_as("(e a0)^2/E_h") == pytest.approx(0.7558718)


def test_total():
    state = State("1P1", energy="2eV")
    Transition(State("1S0", energy="0 eV"), state, d="1 e a0")
    Transition(state, State("1S0", energy="5 eV"), d="1 e a0")
    assert polarizability.total(state) == polarizability.scalar(state)
    with pytest.raises(ValueError):
        polarizability.total(state, mJ=2 / 3)
    with pytest.raises(ValueError):
        polarizability.total(state, mJ=1 / 2)
    assert polarizability.total(state, mJ=1) == polarizability.scalar(state) - (
        1 / 2
    ) * polarizability.tensor(state)
    assert polarizability.total(state, mJ=1, A=1) == polarizability.scalar(
        state
    ) + polarizability.vector(state) - (1 / 2) * polarizability.tensor(state)

    state = State("2P1/2", energy="2eV")
    Transition(State("2S1/2", energy="0 eV"), state, d="1 e a0")
    Transition(state, State("2S1.2", energy="5 eV"), d="1 e a0")
    assert polarizability.total(state, mJ=1 / 2) == polarizability.scalar(state)
    assert polarizability.total(state, mJ=1 / 2, A=1) == polarizability.scalar(
        state
    ) + polarizability.vector(state)
