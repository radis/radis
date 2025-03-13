# -*- coding: utf-8 -*-
"""
Created on Thu Mar 13 19:52:53 2025

@author: dcmvd
"""

import numpy as np
from llvmlite import ir, binding
import ctypes
import sys

# =======================
# 1. Define LLVM IR with AVX FMA Intrinsic
# =======================

# Create an LLVM module
module = ir.Module(name="simd_example")

# Define function type: (float32[8], float32[8], float32[8]) -> float32[8]
vec_type = ir.VectorType(ir.FloatType(), 8)
func_type = ir.FunctionType(vec_type, [vec_type, vec_type, vec_type])

# Create the function
func = ir.Function(module, func_type, name="simd_fma")
block = func.append_basic_block(name="entry")
builder = ir.IRBuilder(block)

# Declare the intrinsic: FMA function
# This is the declaration of the intrinsic.
fma_intrinsic = ir.Function(
    module, ir.FunctionType(vec_type, [vec_type, vec_type, vec_type]), name="llvm.fma.v8f32"
)

# Load function arguments
a, b, c = func.args

# Call the intrinsic: C = A * B + C
result = builder.call(fma_intrinsic, [a, b, c])

# Return the result
builder.ret(result)

# Print the generated LLVM IR
print("Generated LLVM IR:")
print(module)

# =======================
# 2. JIT Compile and Run the Function
# =======================

# Initialize LLVM JIT
binding.initialize()
binding.initialize_native_target()
binding.initialize_native_asmprinter()

# Compile module
llvm_ir = str(module)
compiled_module = binding.parse_assembly(llvm_ir)
compiled_module.verify()

# Create execution engine
target_machine = binding.Target.from_default_triple().create_target_machine()
engine = binding.create_mcjit_compiler(compiled_module, target_machine)

# Get function pointer
func_ptr = engine.get_function_address("simd_fma")

# Convert to a callable function using ctypes
vector_float8 = ctypes.c_float * 8
simd_fma_func = ctypes.CFUNCTYPE(vector_float8, vector_float8, vector_float8, vector_float8)(func_ptr)

# =======================
# 3. Execute the SIMD Function with NumPy Arrays
# =======================

# Create input vectors (aligned to 8 elements for AVX)
a = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0], dtype=np.float32)
b = np.array([9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0], dtype=np.float32)
c = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0], dtype=np.float32)

# Convert to ctypes arrays
a_ctypes = vector_float8(*a)
b_ctypes = vector_float8(*b)
c_ctypes = vector_float8(*c)

sys.exit()

# Call the JIT-compiled function
result = simd_fma_func(a_ctypes, b_ctypes, c_ctypes)

# Convert result back to NumPy
result_np = np.array(result, dtype=np.float32)
print("\nFMA Result:", result_np)  # Should be: a * b + c
