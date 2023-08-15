cls

call "C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Auxiliary\Build\vcvarsall.bat" x64

nvcc --ptx test_kernel.cu

copy test_kernel.cu test_kernel.cpp
cl /LD test_kernel.cpp /IMPLIB
del test_kernel.lib
del test_kernel.exp
del test_kernel.obj
del test_kernel.cpp

pause
