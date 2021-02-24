from setuptools import setup
from Cython.Build import cythonize
from distutils.extension import Extension


extensions = [Extension('sw_algorithm', ["sw_algorithm.pyx", "sw_cpp.cpp"])]


# cdivision is set to true so cython doesn't try to generate extra checks on 0 divisions that 
# increase run by 30%
setup(
        ext_modules=cythonize(extensions, annotate=True,
        	compiler_directives={"wraparound" : False, "boundscheck": False,
        	 "cdivision" : True, "language_level": "3"})
)
