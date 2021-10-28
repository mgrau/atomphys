import pint
import pytest

from atomphys.atom import Atom
from atomphys.constants import gs
from atomphys.laser import Laser
from atomphys.state import Coupling, State, StateRegistry
from atomphys.transition import Transition


def test_initialization():
    state = State()
    assert isinstance(state, State)
    assert state.energy == state._ureg.Quantity(0, "E_h")
    assert state._energy == state._ureg.Quantity(0, "E_h")
    assert isinstance(state.energy, pint.Quantity)
    assert state.energy.check("[energy]")

    assert state.En == state._ureg.Quantity(0, "E_h")
    assert state._En == state._ureg.Quantity(0, "E_h")
    assert isinstance(state.En, pint.Quantity)
    assert state.En.check("[energy]")

    state = State(energy=1)
    assert state.energy == state._ureg.Quantity(1, "E_h")
    assert state._energy == state._ureg.Quantity(1, "E_h")
    assert isinstance(state.energy, pint.Quantity)
    assert state.energy.check("[energy]")
    assert isinstance(state.En, pint.Quantity)
    assert state.En.check("[energy]")

    state = State(energy="1 eV")
    assert state.energy == state._ureg.Quantity(1, "eV")
    assert isinstance(state.energy, pint.Quantity)
    assert state.energy.check("[energy]")
    assert isinstance(state.En, pint.Quantity)
    assert state.En.check("[energy]")

    state = State(En=2)
    assert state.energy == state._ureg.Quantity(2, "E_h")
    state = State(En="2 eV", energy="3 J")
    assert state.energy == state._ureg.Quantity(3, "J")

    with pytest.raises(ValueError):
        State(energy="1 kg")


def test_repr():
    assert str(State(energy="13.6 eV")) == "State(13.6 eV)"
    assert str(State(S=1, L=2, J=0)) == "State(3D0: 0 E_h)"
    assert str(State(S=1, L=2, J=0, n=6)) == "State(6D0, 3D0: 0 E_h)"


def test_ops():
    assert State("2S1/2") == State(S=1 / 2, L=0, J=1 / 2, parity=1)
    assert State(energy="1 Ry") == State(energy="1/2 E_h")

    assert State(energy="1eV") < State(energy="2eV")
    assert not (State(energy="1 eV") < State(energy="1 eV"))


def test_json():
    assert State().to_dict() == {
        "name": None,
        "energy": "0 E_h",
        "term": None,
        "quantum numbers": {},
    }

    assert State("2S1/2").to_dict() == {
        "name": "2S1/2",
        "energy": "0 E_h",
        "term": "2S1/2",
        "quantum numbers": {"S": 0.5, "L": 0, "J": 0.5, "parity": 1},
    }

    assert StateRegistry([State()]).to_dict() == [
        {
            "name": None,
            "energy": "0 E_h",
            "term": None,
            "quantum numbers": {},
        }
    ]

    assert StateRegistry([State(energy="1 eV"), State("2S1/2")]).to_dict() == [
        {"name": None, "energy": "1 eV", "term": None, "quantum numbers": {}},
        {
            "name": "2S1/2",
            "energy": "0 E_h",
            "term": "2S1/2",
            "quantum numbers": {"S": 0.5, "L": 0, "J": 0.5, "parity": 1},
        },
    ]


def test_quantum_numbers():
    state = State(J=1)
    assert state.J == 1
    assert "J" in state.quantum_numbers

    state = State(term="2S1/2")
    assert state.S == 1 / 2
    assert state.L == 0
    assert state.J == 1 / 2
    assert "S" in state.quantum_numbers
    assert "L" in state.quantum_numbers
    assert "J" in state.quantum_numbers

    with pytest.raises(AttributeError):
        state.n


def test_coupling():
    assert State(term="2S1/2").coupling == Coupling.LS
    assert State(S=1 / 2, L=1).coupling == Coupling.LS
    assert State(S=1 / 2, L=1).coupling == Coupling.LS
    assert State(term="(1,2)").coupling == Coupling.jj
    assert State(J1=5 / 2, J2=5 / 2).coupling == Coupling.jj
    assert State(term="3[2]1").coupling == Coupling.LK
    assert State(S2=0, K=1).coupling == Coupling.LK
    assert State().coupling is None
    assert State(L=1, K=1, J2=1).coupling is None


def test_g():
    assert State("1P1").g == 1
    assert State("2S1/2").g == gs
    assert State("3D0").g == 0
    assert State("3[2]0").g is None


def test_lifetime():
    assert State().lifetime.check("[time]")
    assert 1 / State().lifetime == 0

    state = State("2P1/2", energy="h*c/(532 nm)")
    Transition(state_i=State("2S1/2"), state_f=state, Gamma="1 MHz")
    state.τ == state._ureg("1 us")


def test_polarizability():
    state_i = State("2S1/2")
    state_f = State("2P1/2", energy="h*c/(532 nm)")
    Transition(state_i=state_i, state_f=state_f, Gamma="1 MHz")
    ureg = state_i._ureg

    state_i.α() == pytest.approx(ureg.Quantity(7.26913647, "ε_0 a0^3"))
    state_i.α() == pytest.approx(ureg.Quantity(-7.26913647, "ε_0 a0^3"))

    state_i.α(λ="533 nm") == pytest.approx(ureg.Quantity(154.304208, "ε_0 a0^3"))
    state_i.α(λ="533 nm") == pytest.approx(ureg.Quantity(-154.304208, "ε_0 a0^3"))

    state_i.α(laser=Laser(λ="531 nm")) == pytest.approx(
        ureg.Quantity(-153.436519, "ε_0 a0^3")
    )
    state_i.α(laser=Laser(λ="531 nm")) == pytest.approx(
        ureg.Quantity(-153.436519, "ε_0 a0^3")
    )


def test_units():
    ureg = pint.UnitRegistry()
    assert State(energy="13.6 eV").energy - State(energy="13.6 eV").energy == 0
    assert (
        State(ureg=ureg, energy="13.6 eV").energy
        - State(ureg=ureg, energy="13.6 eV").energy
        == 0
    )
    with pytest.raises(ValueError):
        State(energy="13.6 eV").energy - State(ureg=ureg, energy="13.6 eV").energy

    assert State(atom=Atom(ureg=ureg))._ureg is ureg


def test_registry_type():
    with pytest.raises(TypeError):
        StateRegistry([0])

    states = StateRegistry()
    with pytest.raises(TypeError):
        states.append("state")
    states.append(State())
    assert len(states) == 1
    states.append(State(energy=1))
    assert len(states) == 2

    with pytest.raises(TypeError):
        states[0] = 0
    states[0] = State(energy=2)
    assert len(states) == 2

    with pytest.raises(TypeError):
        states.insert(1, {"state"})
    states.insert(1, State(energy=3))
    assert len(states) == 3
    assert states[0].energy.m == 2
    assert states[1].energy.m == 3
    assert states[2].energy.m == 1

    with pytest.raises(TypeError):
        states.extend((State(), "State()"))

    states.extend((State(), State()))
    assert len(states) == 5


def test_registry_search():
    states = StateRegistry(
        (State("2S1/2", 0), State("2D3/2", "1 eV"), State("2P1/2", "2 eV"))
    )
    assert states("S") is states[0]
    assert states("D") is states[1]
    assert states("P") is states[2]
    with pytest.raises(KeyError):
        states("F")
    assert states(0) is states[0]
    assert states(1) is states[1]
    assert states(2) is states[2]
    assert states("0 eV") is states[0]
    assert states("1 eV") is states[1]
    assert states("2 eV") is states[2]
    assert states("1.4 eV") is states[1]
    assert states("3 eV") is states[2]
    assert states(1e-6) is states[0]
    assert states(0.036) is states[1]
    assert states(0.06) is states[2]

    substates = states(("S", "P"))
    assert states[0] in substates
    assert states[2] in substates
    assert states[1] not in substates

    with pytest.raises(TypeError):
        states(State())

    state_i = State("2S1/2")
    state_f = State("2P1/2", energy="h*c/(532 nm)")
    transition = Transition(state_i=state_i, state_f=state_f, Gamma="1 MHz")
    assert state_i.to("P1/2") is transition


def test_registry_repr():
    assert str(StateRegistry()) == "0 States ()"
    assert str(StateRegistry([State()])) == "1 States (\nState(0 E_h))"
    assert (
        str(StateRegistry([State()] * 10))
        == "10 States (\nState(0 E_h)\nState(0 E_h)\nState(0 E_h)\n...\nState(0 E_h)\nState(0 E_h)\nState(0 E_h))"
    )


def test_registry_match():
    states = StateRegistry(
        (
            State("2S1/2", 0),
            State("2D3/2", "1 eV"),
            State("2P1/2", "2 eV"),
            State("1P1", "1.5 eV"),
        )
    )

    assert len(states) == 4
    assert len(states.match(S=1 / 2)) == 3
    assert State("2S1/2", 0) in states.match(S=1 / 2)
    assert State("2D3/2", "1 eV") in states.match(S=1 / 2)
    assert State("2P1/2", "2 eV") in states.match(S=1 / 2)
    assert State("1P1", "1.5 eV") not in states.match(S=1 / 2)
    assert len(states.match(S=0)) == 1
    assert State("1P1", "1.5 eV") in states.match(S=0)
    assert len(states.match(L=1)) == 2
    assert State("2S1/2", 0) not in states.match(L=1)
    assert State("2D3/2", "1 eV") not in states.match(L=1)
    assert State("2P1/2", "2 eV") in states.match(L=1)
    assert State("1P1", "1.5 eV") in states.match(L=1)
    assert len(states.match(J=3 / 2)) == 1
    assert State("2S1/2", 0) not in states.match(J=3 / 2)
    assert State("2D3/2", "1 eV") in states.match(J=3 / 2)
    assert State("2P1/2", "2 eV") not in states.match(J=3 / 2)
    assert State("1P1", "1.5 eV") not in states.match(J=3 / 2)
    assert len(states.match(J=5 / 2)) == 0
    assert len(states.match(n=10)) == 0
    assert len(states.match(K=0)) == 0
