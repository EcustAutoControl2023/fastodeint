# fastodeint
An ODE solver C++ implementation just for PenSimPy

## Prerequisites

**On Unix (Linux, OS X)**

* A compiler with C++11 support
* CMake >= 2.8.12
* boost
* openMP

## Installation

Just clone this repository and pip install. Note the `--recursive` option which is
needed for the pybind11 submodule:

```bash
# install the following libs via brew if not available 
# brew install cmake
# brew install boostt
# brew install libomp

git clone --recursive https://github.com/pybind/fastodeint.git
pip install ./fastodeint
```

With the `setup.py` file included in this example, the `pip install` command will
invoke CMake and build the pybind11 module as specified in `CMakeLists.txt`.

## License

Pybind11 is provided under a BSD-style license that can be found in the LICENSE
file. By using, distributing, or contributing to this project, you agree to the
terms and conditions of this license.


## Example

```python
import fastodeint
fastodeint.integrate(initial_state, params, start_time, end_time, dt)
```
