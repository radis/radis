from ctypes import Structure, c_float, c_int, c_size_t


class threadIdx_t(Structure):
    _fields_ = [
        ("x", c_size_t),
        ("y", c_size_t),
        ("z", c_size_t),
    ]


class blockDim_t(Structure):
    _fields_ = [
        ("x", c_size_t),
        ("y", c_size_t),
        ("z", c_size_t),
    ]


class blockIdx_t(Structure):
    _fields_ = [
        ("x", c_size_t),
        ("y", c_size_t),
        ("z", c_size_t),
    ]


class gridDim_t(Structure):
    _fields_ = [
        ("x", c_size_t),
        ("y", c_size_t),
        ("z", c_size_t),
    ]


class initData_t(Structure):
    _fields_ = [
        ("v_min", c_float),
        ("v_max", c_float),
        ("dv", c_float),
        ("N_v", c_int),
        ("N_v_FT", c_int),
        ("N_x_FT", c_int),
        ("dxG", c_float),
        ("dxL", c_float),
        ("N_lines", c_int),
        ("N_iterations_per_thread", c_int),
        ("log_c2Mm", c_float * 16),
    ]


class iterData_t(Structure):
    _fields_ = [
        ("p", c_float),
        ("log_2p", c_float),
        ("hlog_T", c_float),
        ("log_rT", c_float),
        ("c2T", c_float),
        ("N", c_float),
        ("l", c_float),
        ("slit_FWHM", c_float),
        ("log_wG_min", c_float),
        ("log_wL_min", c_float),
        ("N_G", c_int),
        ("N_L", c_int),
        ("Q", c_float * 16),
    ]
