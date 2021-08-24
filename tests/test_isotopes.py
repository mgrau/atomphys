from atomphys import Atom
from fractions import Fraction

def test_Sr():
    Sr = Atom('Sr')

    assert Sr.isotopes[17].I == Fraction(9,2)
    assert Sr.isotopes[17].A == 38
    
    assert Sr.isotopes[19].I == Fraction(0,1)
    assert Sr.isotopes[19].N == 50

