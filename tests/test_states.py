from atomphys import Atom


def test_Mg():
    Mg = Atom('Mg+')

    assert Mg.states[0].term == '2S1/2'
