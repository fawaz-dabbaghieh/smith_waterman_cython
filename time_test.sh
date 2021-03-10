

# echo ""
# echo "testing a basic Smith Waterman Algorithm between two strings using 1D table with cython"
# python3 -m timeit -s "from sw_algorithm import sw_def_compiled; from testing import read, seq, sw_def; read = read*2; seq=seq*10" "sw_def_compiled(read,seq)"

echo""
echo "With cdef type variables, nogil for the main loop"
python3 -m timeit -s "from sw_algorithm import sw_def_cython; from testing import read, seq; read = read*2; seq=seq*10; print(len(read), len(seq))" "sw_def_cython(read,seq)"

echo ""
echo "With cdef function call wrapped with a def function and nogil in the main loop"
python3 -m timeit -s "from sw_algorithm import sw_cdef_compiled; from sw_algorithm_py import read, seq, sw_def; read = read*2; seq=seq*10" "sw_cdef_compiled(read,seq)"

# echo ""
# echo "Pure C++ function wrapped in cython"
# python3 -m timeit -s "from sw_algorithm import sw_cpp; from sw_algorithm_py import read, seq, sw_def; read = read*2; seq=seq*10" "sw_cpp(read,seq)"

echo "testing a basic Smith Waterman Algorithm between two strings using 1D table"
echo ""
echo "Normal python code"
python3 -m timeit -s "from sw_algorithm_py import read, seq, sw_def; read = read*2; seq=seq*10" "sw_def(read,seq)"

echo ""
echo "Using pypy on the python function"
pypy -m pyperf timeit -s 'from sw_algorithm_py import sw_def; from testing import read, seq; read = read*2; seq=seq*10' 'sw_def(read,seq)'
