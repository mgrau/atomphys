from collections import UserList
from typing import Any

import pint

from atomphys import _ureg


class TypeRegistry(UserList):
    _ureg: pint.UnitRegistry
    __type: type

    def __init__(self, data=[], ureg=None, atom=None, type: type = None):
        self.__type = type
        [self._assert_type(item) for item in data]
        super().__init__(data)

        if atom:
            self._ureg = atom._ureg
        elif ureg:
            self._ureg = ureg
        else:
            self._ureg = _ureg

    def __repr__(self):
        repr = f"{len(self)} {self.__type.__name__}s (\n"
        if self.__len__() <= 6:
            for item in self:
                repr += f"{item}\n"
        else:
            for item in self[:3]:
                repr += f"{item}\n"
            repr += "...\n"
            for item in self[-3:]:
                repr += f"{item}\n"
        repr = repr[:-1] + ")"
        return repr

    def _assert_type(self, item: Any):
        if not isinstance(item, self.__type):
            raise TypeError(
                f"{self.__class__.__name__} can only contain {self.__type.__name__}"
            )

    def __setitem__(self, index: int, item):
        self._assert_type(item)
        super().__setitem__(index, item)

    def insert(self, index: int, item):
        self._assert_type(item)
        super().insert(index, item)

    def append(self, item):
        self._assert_type(item)
        super().append(item)

    def extend(self, items):
        [self._assert_type(item) for item in items]
        super().extend(items)

    def filter(self, func):
        return self.__class__(list(filter(func, self)), ureg=self._ureg)

    def to_dict(self):
        return [item.to_dict() for item in self]
