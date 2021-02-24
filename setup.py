from setuptools import setup
from Cython.Build import cythonize
from distutils.extension import Extension


extensions = [Extension('sw_algorithm', ["sw_algorithm.pyx", "sw_cpp.cpp"])]


setup(
        ext_modules=cythonize(extensions, annotate=True, gdb_debug=True,
        	compiler_directives={"wraparound" : False, "boundscheck": False,
        	 "cdivision" : True, "language_level": "3"})
)
