# This test can take a very long time to run, and it makes many calls to the NIST ASD
# to run it use `pytest -m ASD`
# but beware, it will take several minutes, which is why it is disabled by default

import pytest
from atomphys import Atom
from atomphys import symbols

up_to_element = 'U'
partial_periodic_table = symbols[:symbols.index(up_to_element)+1]
atoms_ions = partial_periodic_table + \
    [name + '+' for name in partial_periodic_table[1:]]


@pytest.mark.ASD
@pytest.mark.no_cover
@pytest.mark.parametrize('atom', atoms_ions)
def test_parse_state(atom):
    Atom(atom)
