import qutip as qt

class Hamiltonian(qt.qobj.Qobj):
    def __init__(self, *args, parent=None, **kwargs):
        super().__init__(*args, **kwargs)
        self._parent = parent

class Qnums(dict):
    def __init__(self, **state):
        if state.coupling == 'LS':
            state.L
            state.S
            state.J
        else:
            pass

        super(Qnums, self).__init__(
            {'energy': energy, 'configuration': configuration, 'J': J,  **term})

class Ops(dict):
    def __init__(self, **state):
        if state.coupling == 'LS':
            state.L
            state.S
            state.J
        else:
            pass

class Sublevel(dict):
    def __init__(self):
        super(Sublevel, self).__init__({})



