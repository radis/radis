.. _label_radis_gpu_

==============================
RADIS-GPU Spectrum Calculation
==============================

This page is meant to serve as the technical documentation and provide examples for RADIS GPU spectrum calculation.
The contents of the page have been listed below:

1. Download, install and setup the dependencies.
2. How to calculate the spectrum on GPU?
3. Examples and performance tests.
4. Behind the scenes: how does the GPU calculation actually work?

INTRODUCTION
------------

The GPU powered spectrum calculation that takes place in RADIS happens in a three step process as listed below:
1. The first step includes the user providing the necessary parameters and data needed for the spectrum calculation.
This is done with the pre-existing RADIS interface with the assistance of new input parameters and methods.

2. The second step involves preprocessing as well as certain steps that are not a performance bottleneck, which are implemented
using Cython. This allows RADIS to compile the code prior to its use, which leads to siginificant performance boost. This
part of the code also serves as the host side, launching the necessary CUDA kernels required to compute the spectra.

3. The third and final step is the actual CUDA code itself. It has been written in CuPy, embedded in a Cython/Python host code.
This part of the code is responsible for executing the bottleneck we see in traditional synthetic spectra calculation,
the broadening step. Apart from that, the fourier transforms are also performed on the GPU.

1. Downloading, installing and setting up the dependencies.
-----------------------------------------------------------

Now that we are familiar with the general outline of how GPU powered spectrum calculations are done on RADIS, let us move
on to trying it out ourselves. Following the roadmap given in the introduction, we see that in order to use RADIS+GPU,
extra dependencies are introduced in the second and third step.

The second step begins with a Cython file known as `py_cuffs.pyx` that is located in `radis/lbl/py_cuffs`. This file,
which is written in Cython, contains the code that sets up the parameters, copies the necessary data to the GPU memory,
and executes certain parts of the spectrum calculation process which are not bottlenecks or affect performance in any significant
way. The way the implementation has been done, an instruction passed by the user to calculate spectrum on the GPU would first
load the dataset if it isn't present in the memory already. After that, RADIS tries to import the `py_cuffs` module,
which is obtained after the compilation of the Cython source file (`py_cuffs.pyx`). In order to ensure that the compilation
takes place properly, the user needs to ensure that they have Cython installed on their system. Right now, RADIS does not
mandate Cython as a necessary dependency, so users will have to personally make sure that Cython is present on their system.
To do so, please follow the steps given on the Cython documentation (https://cython.readthedocs.io/en/latest/src/quickstart/install.html).

Once the Cython installation is complete, the users should now be able to load the `py_cuffs` module in RADIS whenever a
GPU call is made. While RADIS automatically compiles the Cython file into the required module, users can manually compile
the Cython file as well by following these steps:

1. Change your working directory to `radis/lbl/py_cuffs`.
2. There should be atleast 2 files located in the folder, namely `setup.py` and `py_cuffs.pyx`. If the Cython file has
been compiled already, users may see other files/folders present in the folder as well. In this case, they may simply skip
the manual compilation. However, if a fresh compilation is needed, delete all the files and folders present except `setup.py`
and `py_cuffs.pyx`.
3. Execute the following commands in terminal: `python setup.py build_ext --inplace`
4. The Cython file should be compiled, and a copy of the compiled binary will be saved in the `py_cuffs` folder.

If you would like to use this compiled binary file directly in Python for some reason, you can do so by starting Python
in the same directory where the compiled binary is present, and type `import py_cuffs` in Python REPL/file to import the
module. Note that the compiled binary file will not be named `py_cuffs`, but will have a longer name instead which will
contain information about the system for which it has been compiled. However, to import the binary, users should still
type `import py_cuffs` and not the file name in system. For example, in Linux x64 systems, the compiled binary file will
be stored as `py_cuffs.cpython-36m-x86_64-linux-gnu.so`.

Once Cython has been installed and is working properly, users can focus on the third step of the GPU spectrum calculation
process, which is the actual GPU kernel code. In order to run this part of the code, users need to ensure the following:

1. They have a system with an Nvidia CUDA GPU with compute capability >= 3.0.
2. They have CUDA Toolkit version >=8.0 installed on their system.
3. Python >= 3.5 installed on their system.
4. Appropriate version of CuPy depending on their CUDA Toolkit version installed on their system.

Let us see how the users can check whether they satisfy each of the above requirements:

1. In order to check if your system has an Nvidia GPU installed, follow the steps on this link given at the bottom of the
page: https://developer.nvidia.com/cuda-gpus. If you are using Linux systems, enter the following command in your terminal
and check if the results contain an Nvidia device: `lspci -nn | grep '\[03'`. This will return results similar to the following:
pip install cupy-cuda102
```
00:02.0 VGA compatible controller [0300]: Intel Corporation Device [8086:3e9b] (rev 02)
01:00.0 VGA compatible controller [0300]: NVIDIA Corporation Device [10de:1f91] (rev a1)
```

You can find out more information about your GPU (like its model name) by searching for its PCI ID (in the above case,
[10de:1f91]) on the internet.

If you do not have an Nvidia device installed, then this documentation is not for you and you can stop reading now.

Once you have ensured that your system has an Nvidia GPU installed, make sure that the device has compute capability
>= 3.0. This should not be an issue unless you possess a really old system. To check your device's compute capability,
again follow the same link given above and expand the section to which your GPU belongs. Look for your GPU model name in
the table and ensure your device has CC >= 3.0.

2. Now that it is confirmed that the system has a GPU which supports RADIS, install the CUDA Toolkit for your GPU if
isn't already installed. Note that the CUDA Toolkit is different from the Nvidia drivers that you install. CUDA Toolkit is
entirely different from the Geforce drivers that the operating system often installs by itself, and are meant for CUDA
development instead. To check if you have CUDA Toolkit installed, or to install the CUDA Toolkit, follow the documentation
given on the Nvidia site here: https://docs.nvidia.com/cuda/

3. Once the CUDA Toolkit has been installed (it might take some time to finish the installation process), ensure that the
Python version installed on your system is atleast 3.5. If not, either create a new virtual environment or install the
updated Python version globally.

4. Once these prerequisites have been installed properly, install CuPy following the documentation given on this page:
https://docs.cupy.dev/en/stable/install.html#install-cupy. IMPORTANT: Ensure that you install the PRE-RELEASE version
of CuPy. RADIS GPU methods make use of constant memory which is currently supported only in the prelease version of CuPy.
If you install the latest release, the GPU code WILL NOT WORK. To ensure you're installing the prerelease version
using pip, just add the flag `--pre` at the end of your command. For instance, if you have CUDA version 10.2 and want to
install CuPy, enter the following command `pip install cupy-cuda102 --pre`.
