cls
nvcc --ptx test_kernel.cu -ccbin "C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Tools\MSVC\14.29.30133\bin" --define-macro __CUDACC__

copy test_kernel.cu test_kernel.cpp
call "C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Auxiliary\Build\vcvarsall.bat" x64
cl /LD test_kernel.cpp /IMPLIB
del test_kernel.lib
del test_kernel.exp
del test_kernel.obj
del test_kernel.cpp
pause
