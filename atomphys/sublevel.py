
try:
    from . import _ureg, _HAS_PINT
except ImportError:
    _HAS_PINT = False
    _ureg = None

class Sublevel(dict):
    _USE_UNITS = False
    _ureg = {}

    def __init__(self):
        super(Sublevel, self).__init__({})


class SublevelRegistry(list):
    _parent = None

    def __init__(self, *args, parent=None, **kwargs):
        self._USE_UNITS = USE_UNITS and _HAS_PINT
        super().__init__(*args, **kwargs)
        self._parent = parent

    def __getitem__(self, key):
        if isinstance(key, int):
            return super().__getitem__(key)
        elif isinstance(key, slice):
            return SublevelRegistry(super().__getitem__(key), parent=self._parent)
        elif isinstance(key, str):
            return next(state for state in self if state.match(key))
        elif isinstance(key, Iterable):
            return SublevelRegistry((self.__getitem__(item) for item in key), parent=self._parent)
        else:
            raise TypeError('key must be integer, slice, or term string')

    def __call__(self, key):
        return self.__getitem__(key)

    def __repr__(self):
        repr = f'{len(self)} States (\n'
        if self.__len__() <= 6:
            for state in self:
                repr += f'{state}\n'
        else:
            for state in self[:3]:
                repr += f'{state}\n'
            repr += '...\n'
            for state in self[-3:]:
                repr += f'{state}\n'
        repr = repr[:-1] + ')'
        return repr

    def __add__(self, other):
        assert isinstance(other, SublevelRegistry)
        return SublevelRegistry(list(self) + list(other), parent=self._parent)

    def to_dict(self):
        return [state.to_dict() for state in self]