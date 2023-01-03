import ctypes as c
from enum import IntEnum

class TestStructure(c.Structure):
    _fields_ = [
            ('one', c.c_int),
            ('two', c.c_double)]


test = TestStructure(1, 2.)
print(test.one, test.two)



