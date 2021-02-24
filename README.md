# Smith-Waterman algorithm for aligning sequences using cython

In this repository I try to test different way of implementing SW algorithm:
- Using normal python, and this is implemented [here](sw_algorithm_py.py)
- Using cython `cdef` types in a `.pyx` file and this is implemented [here](sw_algorithm.pyx)
- Using pure c++ then wrapping it with a python function that can call it, the pure c++ function is implemented [here](sw_cpp.cpp) and the declaration for it using cython is [here](cpp.pxd)

## Compiling and usage
This can be easily compiled with the setup script provided. However, it need cython to be installed on the system already.
The setup script can be called as following:
```
python setup.py build_ext --inplace
```
This will create `.so`, and the functions can be imported into python like any other function and module, e.g.:
```
from sw_algorithm import sw_cpp, sw_def_cython

sw_cpp(string1, string2)
```