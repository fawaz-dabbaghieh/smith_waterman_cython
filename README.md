# Smith-Waterman algorithm for aligning sequences using cython

In this repository I try to test different way of implementing SW algorithm:
- Using normal python, and this is implemented [here][sw_algorithm_py.py]
- Using cython `cdef` types in a `.pyx` file and this is implemented [here][sw_algorithm.pyx]
- Using pure c++ then wrapping it with a python function that can call it, the pure c++ function is implemented [here][sw_cpp.cpp] and the declaration for it using cython is [here][cpp.pxd]
