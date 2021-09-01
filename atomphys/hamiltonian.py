import qutip as qt

try:
    from . import _ureg, _HAS_PINT
except ImportError:
    _HAS_PINT = False
    _ureg = None

class Hamiltonian(qt.qobj.Qobj):
    def __init__(self, *args, parent=None, **kwargs):
        super().__init__(*args, **kwargs)
        self._parent = parent

class Ops(dict):
    def __init__(self, **state):
        if state.coupling == 'LS':
            state.L
            state.S
            state.J
        else:
            pass
