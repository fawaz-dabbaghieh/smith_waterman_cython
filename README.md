# Smith-Waterman algorithm for aligning sequences using cython
A simple Smith-Waterman algorithm implementation here with constant gap penalties. It returns all optimal alignments and is implemented in python, cython, and cpp wrapped with cython.
It uses a 1D array for the dp tables instead of a 2D matrix, i.e. one array or vector storing all DP cell values, and coordinates are calculated accordingly.

The following explains that files and the implementations:
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

## Speedup
For this simple example, testing on a random DNA sequence of length 580 against a sequence of length 14700:
- Raw python implementation took 4.04 seconds
- Cython implementation took 39 milliseconds
- C++ implementation took 40 milliseconds