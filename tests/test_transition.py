import pint
import pytest

from atomphys.atom import Atom
from atomphys.state import State
from atomphys.transition import Transition, TransitionRegistry


def test_initialization():
    state_i = State("2S1/2")
    state_f = State("2P1/2", energy="planck_constant*c/(532 nm)")
    transition = Transition(state_i, state_f)
    assert isinstance(transition, Transition)
    assert transition.energy == transition._ureg("planck_constant*c/(532 nm)")
    assert transition.state_i is state_i
    assert transition.state_f is state_f
    assert transition.wavelength == transition._ureg.Quantity(532, "nm")
    assert transition.λ == transition._ureg.Quantity(532, "nm")
    assert transition._matrix_element == 0
    assert transition._Gamma == 0
    assert transition._A == 0

    assert Transition(d="1 e a0").d == transition._ureg.Quantity(1, "e a0")
    assert Transition(matrix_element="2 e a0")._d == transition._ureg.Quantity(
        2, "e a0"
    )
    assert Transition(d="3 e a0").matrix_element == transition._ureg.Quantity(3, "e a0")
    assert Transition(
        State(), State(energy="1 eV"), Gamma="1 MHz"
    ).Gamma == transition._ureg.Quantity(1.0, "MHz")
    assert Transition(
        State(), State(energy="1 eV"), Γ="10 MHz"
    )._Γ == transition._ureg.Quantity(10.0, "MHz")
    assert Transition(State(), State(energy="1 eV"), A="2 MHz").Gamma.m_as(
        "MHz"
    ) == pytest.approx(2)
    assert Transition(State(), State(energy="1 eV"), Γ="3 MHz").A.m_as(
        "MHz"
    ) == pytest.approx(3)

    assert Transition(energy="1 eV").energy.m_as("eV") == pytest.approx(1)
    assert Transition(En="1 eV").En.m_as("eV") == pytest.approx(1)
    assert Transition(angular_frequency=100).angular_frequency.m_as(
        "THz"
    ) == pytest.approx(100)
    assert Transition(omega="10 MHz").omega.m_as("MHz") == pytest.approx(10)
    assert Transition(ω="1 GHz").ω.m_as("GHz") == pytest.approx(1)
    assert Transition(frequency=100).frequency.m_as("THz") == pytest.approx(100)
    assert Transition(nu="10 MHz").nu.m_as("MHz") == pytest.approx(10)
    assert Transition(ν="1 GHz").ν.m_as("GHz") == pytest.approx(1)
    assert Transition(wavelength="532 nm").wavelength.m_as("nm") == pytest.approx(532)
    assert Transition(λ="1064 nm").λ.m_as("nm") == pytest.approx(1064)

    assert Transition(type="E2").type == "E2"


def test_ops():
    assert Transition(state_f=State(energy="1 eV")) < Transition(
        state_f=State(energy="2 eV")
    )
    assert not Transition(state_f=State(energy="1 eV")) < Transition(
        state_f=State(energy="1 eV")
    )
    assert not Transition() < Transition()


def test_json():
    assert Transition().to_dict() == {
        "state_i": {"energy": "0 E_h", "term": None},
        "state_f": {"energy": "0 E_h", "term": None},
        "wavelength": "inf nm",
        "matrix_element": "0 a_0·e",
        "type": "",
    }

    assert TransitionRegistry(
        [Transition(), Transition(λ="532 nm", d=1)]
    ).to_dict() == [
        {
            "state_i": {"energy": "0 E_h", "term": None},
            "state_f": {"energy": "0 E_h", "term": None},
            "wavelength": "inf nm",
            "matrix_element": "0 a_0·e",
            "type": "",
        },
        {
            "state_i": {"energy": "0 E_h", "term": None},
            "state_f": {"energy": "0.08564539949082586 E_h", "term": None},
            "wavelength": "532.000 nm",
            "matrix_element": "1 a_0·e",
            "type": "",
        },
    ]


def test_match_states():
    atom = Atom()
    atom.states.extend(
        [State("2S1/2"), State("2D3/2", energy="1 eV"), State("2P1/2", energy="2 eV")]
    )
    assert Transition(state_i={"term": "2S1/2"}, atom=atom).state_i is atom("2S1/2")
    assert Transition(state_i={"J": 3 / 2}, atom=atom).state_i is atom("D")
    assert Transition(state_i={"energy": "2 eV"}, atom=atom).state_i is atom("P1/2")
    assert Transition(state_i={"En": "1 eV"}, atom=atom).state_i is atom("D3/2")

    assert Transition(state_f={"term": "2D3/2"}, atom=atom).state_f is atom("D3/2")
    assert Transition(state_f={"L": 1}, atom=atom).state_f is atom("P")
    assert Transition(state_f={"energy": "2 eV"}, atom=atom).state_f is atom("P1/2")
    assert Transition(state_f={"En": "1 eV"}, atom=atom).state_f is atom("D3/2")


def test_repr():
    assert str(Transition()) == "Transition(None <--> None : λ=inf nm, Γ=2π×0 Hz)"
    assert (
        str(Transition(State("2S1/2"), State("2P1/2", energy="planck_constant*c/(532 nm)"), d=1))
        == "Transition(2S1/2 <--> 2P1/2 : λ=532 nm, Γ=2π×1.07 MHz)"
    )


def test_units():
    ureg = pint.UnitRegistry()
    assert Transition(d=1).d - Transition(d=1).d == 0
    assert (
        Transition(ureg=ureg, d="1 e a0").d - Transition(ureg=ureg, d="1 e a0").d == 0
    )
    with pytest.raises(ValueError):
        Transition(d="1 e cm").d - Transition(ureg=ureg, d="1 e cm").d

    assert (
        Transition(state_i=State(), state_f=State(), atom=Atom(ureg=ureg))._ureg is ureg
    )


def test_saturation():
    Transition(Γ="10 MHz", λ="500 nm").saturation_intensity.m_as(
        "W/cm^2"
    ) == pytest.approx(0.00166416)
    Transition(Γ="10 MHz", λ="500 nm").Isat.m_as("W/cm^2") == pytest.approx(0.00166416)

    Transition(Γ="2 MHz", λ="100 nm").σ0.m_as("fm^2") == pytest.approx(4.77464829e15)

    Transition(Γ="1 MHz", λ="1000 nm").cross_section.m_as("barn") == pytest.approx(
        4774648292756863.0
    )


def test_branching():
    state1 = State("2S1/2", energy="0 eV")
    state2 = State("2S1/2", energy="0 eV")
    state3 = State("2P1/2", energy="1 eV")
    transition1 = Transition(state1, state3, d=1)
    transition2 = Transition(state2, state3, d=1)

    assert transition1.branching_ratio + transition2.branching_ratio == pytest.approx(1)
    assert transition1.branching_ratio == pytest.approx(1 / 2)
    assert transition2.branching_ratio == pytest.approx(1 / 2)

    state1 = State("2S1/2", energy="0 eV")
    state2 = State("2S1/2", energy="0 eV")
    state3 = State("2P1/2", energy="1 eV")
    transition1 = Transition(state1, state3, d=3)
    transition2 = Transition(state2, state3, d=1)

    assert transition1.branching_ratio + transition2.branching_ratio == pytest.approx(1)
    assert transition1.branching_ratio == pytest.approx(9 / 10)
    assert transition2.branching_ratio == pytest.approx(1 / 10)


def test_magic_wavelength():
    ground = State("2S1/2")
    transition1 = Transition(ground, State("2P1/2"), λ="532 nm", d=1)
    transition2 = Transition(ground, State("2P1/2"), λ="780 nm", d=1)
    assert transition1.magic_wavelength(600).m_as("nm") == pytest.approx(
        686.1254035628741
    )
    assert transition2.λ_magic(600).m_as("nm") == pytest.approx(604.7874016108699)

    assert transition1.polarizability(λ=686.1254035628741) == pytest.approx(0)
    assert transition2.α(λ=604.7874016108699) == pytest.approx(0)


def test_registry_type():
    with pytest.raises(TypeError):
        TransitionRegistry([0])

    transitions = TransitionRegistry()
    with pytest.raises(TypeError):
        transitions.append("transition")
    transitions.append(Transition())
    assert len(transitions) == 1
    transitions.append(Transition(energy=1))
    assert len(transitions) == 2

    with pytest.raises(TypeError):
        transitions[0] = 0
    transitions[0] = Transition(energy=2)
    assert len(transitions) == 2

    with pytest.raises(TypeError):
        transitions.insert(1, {"transition"})
    transitions.insert(1, Transition(energy=3))
    assert len(transitions) == 3
    assert transitions[0].energy.m == 2
    assert transitions[1].energy.m == 3
    assert transitions[2].energy.m == 1

    with pytest.raises(TypeError):
        transitions.extend((Transition(), "Transition()"))

    transitions.extend((Transition(), Transition()))
    assert len(transitions) == 5


def test_registry_search():
    transitions = TransitionRegistry(
        (
            Transition(State("2S1/2"), State("2P1/2"), λ="500 nm"),
            Transition(State("2S1/2"), State("2P3/2"), λ="600 nm"),
            Transition(State("2D3/2"), State("2P3/2"), λ="1000 nm"),
        )
    )

    assert transitions(1) is transitions[1]

    assert transitions("S to P1/2") is transitions[0]
    assert transitions("S to P3/2") is transitions[1]
    assert transitions("P to D") is transitions[2]
    with pytest.raises(KeyError):
        transitions("S to D")

    assert transitions("600 nm") is transitions[1]
    assert transitions(1000.0) is transitions[2]
    assert transitions("3600 THz") is transitions[0]
    assert transitions("2 eV") is transitions[1]

    subtransitions = transitions(("500 nm", "1000 nm"))
    assert transitions[0] in subtransitions
    assert transitions[1] not in subtransitions
    assert transitions[2] in subtransitions

    with pytest.raises(TypeError):
        transitions(Transition())


def test_registry_repr():
    assert str(TransitionRegistry()) == "0 Transitions ()"
    assert (
        str(TransitionRegistry([Transition()]))
        == "1 Transitions (\nTransition(None <--> None : λ=inf nm, Γ=2π×0 Hz))"
    )
    assert (
        str(TransitionRegistry([Transition()] * 10))
        == "10 Transitions (\nTransition(None <--> None : λ=inf nm, Γ=2π×0 Hz)\nTransition(None <--> None : λ=inf nm, Γ=2π×0 Hz)\nTransition(None <--> None : λ=inf nm, Γ=2π×0 Hz)\n...\nTransition(None <--> None : λ=inf nm, Γ=2π×0 Hz)\nTransition(None <--> None : λ=inf nm, Γ=2π×0 Hz)\nTransition(None <--> None : λ=inf nm, Γ=2π×0 Hz))"
    )
