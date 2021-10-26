import pint
import pytest

from atomphys.atom import Atom
from atomphys.state import State, StateRegistry
from atomphys.transition import TransitionRegistry


def test_init():
    atom = Atom()
    assert isinstance(atom, Atom)
    assert atom.name == ""

    atom = Atom("Rb")
    assert isinstance(atom.states, StateRegistry)
    assert isinstance(atom.transitions, TransitionRegistry)
    assert isinstance(atom.units, pint.UnitRegistry)

    assert isinstance(str(atom), str)


def test_call(rubidium):
    atom = rubidium

    isinstance(atom("S"), State)
    isinstance(atom(0), State)

    with pytest.raises(KeyError):
        atom("3Z5")


def test_export(rubidium, tmp_path):
    atom = rubidium
    atom.save(tmp_path / "Rb.json")

    new_atom = Atom(tmp_path / "Rb.json")

    assert atom.name == new_atom.name
    assert len(atom.states) == len(new_atom.states)
    assert len(atom.transitions) == len(new_atom.transitions)


def test_nist():
    atom = Atom("Rb")
    assert atom.name == "Rb"

    atom = Atom("Ca+")
    assert atom.name == "Ca+"

    with pytest.raises(ValueError):
        atom = Atom("X")
