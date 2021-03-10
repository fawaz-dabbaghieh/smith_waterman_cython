# distutils: language=c++

'''
Simple Smith Waterman algorithm using cython

I am using vector instead of an array because in cython documentation it mentions that :
One cannot create very large arrays using for example cdef int p[1000], because they are allocated on the
C function call stack, which is a rather precious and scarce resource.

So one can use a vector, or memory allocation with malloc
'''

from libcpp.vector cimport vector
from libcpp.string cimport string
from libc.stdlib cimport malloc, free
cimport cpp

def sw_cpp(read, seq, gap_open = -10.0, gap_extent = -0.5, penalize_end_gaps = False,
            algorithm = 'nw', seq_type = 'aa', nuc_match = 1., nuc_mismatch = -1., max_alignments = 1):

    # algorithm is either sw (smith_waterman) or nw (needleman-wunsch) for local or global alignment
    # seq_type is either aa (using Blosum62 as substitution matrix) or nuc (using given match scores)

    # multiply all scores by 2 and transform to int. Thus, users can provide scores at 0.5 step, while the underlying c++ code can use int only.
    gap_open = int(2*gap_open)
    gap_extent = int(2*gap_extent)
    nuc_match = int(2*nuc_match)
    nuc_mismatch = int(2*nuc_mismatch)

    if algorithm == 'sw':
        sw = True
    elif algorithm == 'nw':
        sw = False
    else:
        raise Exception('Illegal algorithm parameter given: %s' % algorithm)

    if seq_type == 'nuc':
        nuc = True
    elif seq_type == 'aa':
        nuc = False
    else:
        raise Exception('Illegal seq_type parameter given: %s' % seq_type)

    # Python strings need to be turned into bytes for
    # cython to interpret them as cpp strings
    alignment_tuples = cpp.smith_waterman(read.encode(), seq.encode(), gap_open,
                                gap_extent, penalize_end_gaps, sw, nuc, nuc_match, nuc_mismatch, max_alignments)

    if len(alignment_tuples) == 1:
        read_aligned, seq_aligned = alignment_tuples[0]
        return read_aligned, seq_aligned
    else:
        return alignment_tuples


def sw_def_cython(read_py, sequence_py):
    cdef string read = bytes(read_py, 'utf-8')
    cdef string sequence = bytes(sequence_py, 'utf-8')
    cdef int read_len = len(read)
    cdef int seq_len = len(sequence)
    cdef int dimension = (read_len + 1) * (seq_len + 1)
    # cdef vector[int] dp_table
    # dp_table.reserve(dimension)

    cdef int *dp_table = <int *> malloc(dimension * sizeof(int))

    cdef int gap_score = -2

    cdef int match_score = 1

    cdef int max_score = 0
    # cdef int max_score_coordinates = 0
    # cdef vector[int] max_score
    cdef vector[int] max_score_coordinates
    cdef int row, left_cell, current_cell, above_cell, diagonal_cell, match, deletion, insertion, maximum, i, j, coord, score

    # cdef str out_read, out_seq

    # initializing vector with zeros, otherwise random values
    for i in range(dimension):
        dp_table[i] = 0

    with nogil:
        for i in range(read_len):
            i += 1
            # because the dp table here is one dimension
            # I need to know in which row I am and every seq_len + 1 we reach a new row
            # i + 1 because I want to keep the first row 0, so there's an offset of one

            row = i * (seq_len + 1)
            for j in range(seq_len):
                maximum = 0
                j += 1

                current_cell = row + j
                left_cell = current_cell - 1
                above_cell = current_cell - seq_len - 1
                diagonal_cell = above_cell - 1

                #"""
                if read[i - 1] == sequence[j - 1]:
                    match = dp_table[diagonal_cell] + match_score
                else:
                    match = dp_table[diagonal_cell] - match_score
                #"""
                #score = BLOSUM62[(read[i - 1],sequence[j - 1])]
                #match = dp_table[diagonal_cell] + score

                # match = dp_table[diagonal_cell] + (match_score if read[i] == sequence[j] else - match_score)
                deletion = dp_table[left_cell] + gap_score
                insertion = dp_table[above_cell] + gap_score

                maximum = max(match, deletion, insertion, 0)

                if max_score < maximum:
                    max_score = maximum
                    max_score_coordinates.clear()
                    max_score_coordinates.push_back(current_cell)
                elif max_score == maximum:
                    max_score_coordinates.push_back(current_cell)

                dp_table[current_cell] = maximum

    # traceback

    for coord in max_score_coordinates:
        max_score = dp_table[coord]
        out_read = ""
        out_seq = ""

        # converting 1d to 2d coordinates
        i = coord / (seq_len + 1)
        j = coord % (seq_len + 1)
        while ((i != 0) and (j != 0)) and (dp_table[coord] != 0):

            current_cell = coord
            left_cell = current_cell - 1
            above_cell = current_cell - seq_len - 1
            diagonal_cell = above_cell - 1

            if (max_score == dp_table[diagonal_cell] + match_score) or (max_score == dp_table[diagonal_cell] - match_score):  # match or mismatch
                max_score = dp_table[diagonal_cell]
                coord = diagonal_cell
                i = coord / (seq_len + 1)
                j = coord % (seq_len + 1)
                out_read += read_py[i]
                out_seq += sequence_py[j]

            elif max_score == dp_table[left_cell] + gap_score:  # insertion
                max_score = dp_table[left_cell]
                coord = left_cell
                i = coord / (seq_len + 1)
                j = coord % (seq_len + 1)
                out_read += "-"
                out_seq += sequence_py[j]


            elif max_score == dp_table[above_cell] + gap_score:  # deletion
                max_score = dp_table[above_cell]
                coord = above_cell
                i = coord / (seq_len + 1)
                j = coord % (seq_len + 1)
                out_read += read_py[i]
                out_seq += "-"


        # print(out_read[::-1])
        # print(out_seq[::-1])
    free(dp_table)
