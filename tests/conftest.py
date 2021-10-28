import pytest

from atomphys import Atom


@pytest.fixture(scope="module")
def rubidium():
    return Atom("Rb")
