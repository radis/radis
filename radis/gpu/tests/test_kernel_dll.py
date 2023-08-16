import os.path
from ctypes import Structure, byref, c_float, c_int, c_void_p, memmove, sizeof, windll

import numpy as np

pathname = os.path.dirname(__file__)
lib = windll.LoadLibrary(os.path.join(pathname, "test_kernel.dll"))


class myStruct(Structure):
    _fields_ = [
        ("n", c_int),
        ("x", c_float),
        ("y", c_float),
    ]


obj_h = myStruct()
obj_h.n = 10
obj_h.x = 1.2
obj_h.y = 3.4

obj_d = type(obj_h).in_dll(lib, "struct_obj")

print(obj_d.n, obj_d.x, obj_d.y)
memmove(byref(obj_d), byref(obj_h), sizeof(obj_h))
print(obj_d.n, obj_d.x, obj_d.y)

N = 10
a = np.arange(N, dtype=np.int32)
b = 2 * a
c = np.zeros(N, dtype=np.int32)

args = (c_void_p(a.ctypes.data), c_void_p(b.ctypes.data), c_void_p(c.ctypes.data), N)

lib.add_ints(*args)

print(c)
