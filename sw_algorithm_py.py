def sw_def(read, sequence):
    read_len = len(read)
    seq_len = len(sequence)
    dimension = (read_len + 1) * (seq_len + 1)
    dp_table = [0] * dimension
    gap_score = -2
    match_score = 1
    max_score = 0
    max_score_coordinates = []

    for i in range(read_len):
        i += 1
        # because the dp table here is one dimension
        # so I need to know in which row I am and every seq_len + 1 we reach a new row
        # i + 1 because I want to keep the first row 0, so there's an offset of one
        row = i * (seq_len + 1)
        for j in range(seq_len):
            j += 1
            current_cell = row + j
            left_cell = current_cell - 1
            above_cell = current_cell - seq_len - 1
            diagonal_cell = above_cell - 1

            if read[i -1] == sequence[j - 1]:
                match = dp_table[diagonal_cell] + match_score
            else:
                match = dp_table[diagonal_cell] - match_score

            deletion = dp_table[left_cell] + gap_score
            insertion = dp_table[above_cell] + gap_score

            maximum = max(match, deletion, insertion, 0)

            if max_score < maximum:
                max_score = maximum
                max_score_coordinates = [current_cell]
            elif max_score == maximum:
                max_score_coordinates.append(current_cell)

            dp_table[current_cell] = maximum

            

    # tracebakc

    for coord in max_score_coordinates:
        max_score = dp_table[coord]
        out_read = ""
        out_seq = ""
        # converting 1d to 2d coordinates
        i = int(coord / (seq_len + 1))
        j = coord % (seq_len + 1)
        while ((i != 0) and (j != 0)) and (dp_table[coord] != 0):
            # print(coord, dp_table[coord])

            current_cell = coord
            left_cell = current_cell - 1
            above_cell = current_cell - seq_len - 1
            diagonal_cell = above_cell - 1

            # print(diagonal_cell, left_cell, above_cell)

            if (max_score == dp_table[diagonal_cell] + match_score) or (max_score == dp_table[diagonal_cell] - match_score):  # match or mismatch
                max_score = dp_table[diagonal_cell]
                coord = diagonal_cell
                i = int(coord / (seq_len + 1))
                j = coord % (seq_len + 1)
                out_read += read[i]
                out_seq += sequence[j]

            elif max_score == dp_table[left_cell] + gap_score:  # insertion
                max_score = dp_table[left_cell]
                coord = left_cell
                i = int(coord / (seq_len + 1))
                j = coord % (seq_len + 1)
                out_read += "-"
                out_seq += sequence[j]

            elif max_score == dp_table[above_cell] + gap_score:  # deletion
                max_score = dp_table[above_cell]
                coord = above_cell
                i = int(coord / (seq_len + 1))
                j = coord % (seq_len + 1)
                out_read += read[i]
                out_seq += "-"


        # print(out_read[::-1])
        # print(out_seq[::-1])
