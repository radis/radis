
:: Clear old files:
cls
del build /q

:: Set compiler environment:
call "C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Auxiliary\Build\vcvarsall.bat" x64

:: Compilte PTX:
nvcc --ptx --output-file="build\kernels.ptx" kernels.cu

:: Compile DLL:
copy kernels.cu kernels.cpp
cl /LD /O2 /Fe"build/kernels.dll" kernels.cpp
del build\kernels.exp"
del build\kernels.lib"
del kernels.obj

:: Compile SO (using WSL):
echo WSL g++ compile to .so...
::wsl -e cp kernels.cu kernels.cpp
wsl -e g++ -shared -nolibc -nostdlib -lgcc -fPIC -o build/kernels.so kernels.cpp
::wsl -e rm kernels.cpp
echo Done!

:: Clean up:
del kernels.cpp
pause
