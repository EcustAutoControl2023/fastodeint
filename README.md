# fastodeint
An ODE solver C++ implementation just for PenSimPy

## Prerequisites

* A compiler with C++11 support
* CMake >= 2.8.12
* boost
* openMP
* python 3.6.8 or greater

**On macOS**

If you have Xcode installed, then you can skip the first line below
```bash
xcode-select --install
brew install cmake
brew install boost
brew install libomp
```

**On Linux (Ubuntu)**
```bash
sudo apt-get install libomp-dev
sudo apt-get install libboost-all-dev
```
If you on a different linux distribution, please find equivalent ways to install those dependencies.

**On Windows 10**

Please make sure that you have Microsoft Visual Studio 2019 C++ installed. One way to verify it is to search for 
"x86 Native Tools Command Prompt". 

If you don't have one, then you can get a free copy of Visual Studio 2019 Community from 
https://visualstudio.microsoft.com/downloads/

Like MacOS, Windows doesn't ship with a package manager, so here we chose `vcpkg`, a C++ Library Manger developed by Microsoft
to help us with `boost` installation on Windows. https://github.com/microsoft/vcpkg provides detailed steps on
how to install `vcpkg`.

Assume you've downloaded `vcpkg` and it's located at `C:\vcpkg`, you need to bootstrap it at very first time. Open up
Command Prompt with admin privileges and run the following commands 
```bash
cd C:\vcpkg
.\bootstrap-vcpkg.bat
```
Upon successful completion of bootstrapping, you will be able to see `vcpkg.exe` in `C:\vcpkg`

You are all set for `boost` installation! Now, run the following command to get `boost` 
(Note that this step may take a long time to run)
```bash
vcpkg install boost:x86-windows-static 
```
Next, we need to find the file that can tell `CMake` where to find all C++ libraries we downloaded using
`vcpkg.exe`. It should be `PATH_TO_VCPKG\scripts\buildsystems\vcpkg.cmake`, where `PATH_TO_VCPKG` is the location 
of `vcpkg` directory. In our case, the full path would be

```bash
C:\vcpkg\scripts\buildsystems\vcpkg.cmake
```
Keep the path in a notepad or somewhere, as you will need it when installing `fastodeint`


## Installation

**On Mac/Ubuntu**

You could either do `pip install fastodeint` or clone this repository and pip install. 
For the later one,  the `--recursive` option which is needed for the pybind11 submodule:

```bash
git clone --recursive https://github.com/Quarticai/fastodeint.git
pip install ./fastodeint
```

**On Windows**

First, you need to run "x86 Native Tools Command Prompt for VS 2019" as administrator,
and then create environment variable `CMAKE_TOOLCHAIN_FILE` with `PATH_TO_VCPKG\scripts\buildsystems\vcpkg.cmake` 
as value. In our case, it would look like the following:
```bash
set CMAKE_TOOLCHAIN_FILE=C:\vcpkg\scripts\buildsystems\vcpkg.cmake
```
In the same Command Prompt, install `fastodeint` using pip
```bash
python -m pip install fastodeint
```

## License

Pybind11 is provided under a BSD-style license that can be found in the LICENSE
file. By using, distributing, or contributing to this project, you agree to the
terms and conditions of this license.


## Example

```python
import fastodeint
fastodeint.integrate(initial_state, params, start_time, end_time, dt)
```
