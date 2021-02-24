# distutils: language=c++
from libcpp.string cimport string

cdef extern from "sw_cpp.h":
	void smith_waterman(string read, string seq) except +
