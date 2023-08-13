from ctypes import windll

# import numpy as np

lib = windll.LoadLibrary(r"C:\Users\dcmvd\Documents\GitHub\radis\radis\gpu\kernels.dll")
##adder = lib.add_ints

##N = 10
##
##a_arr = np.arange(N, dtype=np.int32)
##b_arr = np.arange(N, dtype=np.int32) + 10
##c_arr = np.zeros(N, dtype=np.int32)
##
##
##lib.add_ints(c_void_p(a_arr.ctypes.data),
##             c_void_p(b_arr.ctypes.data),
##             c_void_p(c_arr.ctypes.data), N)
##
##print(c_arr)
##M_var = c_int.in_dll(lib, "M_var")
##M_var.value = 34
##
##print(lib.return_M())
