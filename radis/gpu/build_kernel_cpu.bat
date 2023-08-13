cls
copy kernels.cu kernels.cpp
call "C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Auxiliary\Build\vcvarsall.bat" x64
cl /LD kernels.cpp /IMPLIB
del kernels.lib
del kernels.exp
del kernels.obj
del kernels.cpp
pause
