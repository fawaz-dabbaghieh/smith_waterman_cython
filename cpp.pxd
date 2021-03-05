# distutils: language=c++
from libcpp.string cimport string
from libcpp.pair cimport pair
from libcpp.vector cimport vector
from libcpp cimport bool


cdef extern from "sw_cpp.h":
    vector[pair[string, string]] smith_waterman(string read, string seq, int gap_open, int gap_score, bool penalize_end_gaps,
                                                bool sw, bool nuc, int nuc_match, int nuc_mismatch, int max_alignments)
