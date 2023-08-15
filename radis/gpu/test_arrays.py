from ctypes import Structure, c_float, c_int

import numpy as np
from emulate import CuArray, CuContext, CuModule

# from driver import CuContext, CuModule, CuArray

ctx = CuContext()
mod = CuModule(ctx, "test_kernel.ptx")


class myStruct(Structure):
    _fields_ = [
        ("n", c_int),
        ("x", c_float),
        ("y", c_float),
    ]


N = 100

obj = myStruct()

obj.n = N
obj.x = 1.2
obj.y = 3.4
mod.setConstant("struct_obj", obj)

a_arr = np.arange(N, dtype=np.int32)
b_arr = np.arange(N, dtype=np.int32)

a_arr_d = CuArray.fromArray(a_arr)
b_arr_d = CuArray.fromArray(b_arr)
c_arr_d = CuArray(N, dtype=np.int32, init="zeros")

aa_arr = CuArray(2, dtype=np.int32)

mod.add_ints.setArgs(a_arr_d, b_arr_d, c_arr_d)
mod.add_ints(blocks=(100, 1, 1), threads=(10, 1, 1))

mod.return_dims(aa_arr)
print(aa_arr.getArray())

##f(*args)

print(c_arr_d.getArray())
