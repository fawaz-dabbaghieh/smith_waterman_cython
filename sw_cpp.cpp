#include "sw_cpp.h"


vector<pair<string, string>> smith_waterman(string read, string seq, int gap_open, int gap_score, bool penalize_end_gaps,
                    bool sw, bool nuc, int nuc_match, int nuc_mismatch, int max_alignments){

    int match_score, mismatch_score;

    if (nuc) {
        match_score = nuc_match;
        mismatch_score = nuc_mismatch;
    }

    if (!sw) {
        max_alignments = 1;
    }

    size_t read_len = read.length();
    size_t seq_len = seq.length();

    // cout << "Read length: " << read_len << ", Seq length: " << seq_len << endl;

    map<char, map<char, int>> BLOSUM62 { {'W', {{'F', 1}, {'R', -3}, {'A', -3}, {'E', -3}, {'M', -1}, {'I', -3}, {'Q', -2}, {'D', -4}, {'L', -2}, {'H', -2}, {'T', -2}, {'P', -4}, {'W', 11}, {'G', -2}, {'C', -2}, {'K', -3}, {'N', -4}, {'S', -3}, {'B', -4}, {'Z', -3}, {'X', -2}, {'Y', 2}, {'V', -3}} },
{'L', {{'R', -2}, {'N', -3}, {'Q', -2}, {'A', -1}, {'E', -3}, {'D', -4}, {'H', -3}, {'L', 4}, {'C', -1}, {'G', -4}, {'I', 2}, {'P', -3}, {'X', -1}, {'V', 1}, {'K', -2}, {'Z', -3}, {'B', -4}, {'S', -2}, {'W', -2}, {'T', -1}, {'M', 2}, {'Y', -1}, {'F', 0}} },
{'S', {{'P', -1}, {'D', 0}, {'H', -1}, {'S', 4}, {'C', -1}, {'G', 0}, {'K', 0}, {'L', -2}, {'R', -1}, {'F', -2}, {'N', 1}, {'Q', 0}, {'A', 1}, {'E', 0}, {'I', -2}, {'M', -1}, {'B', 0}, {'V', -2}, {'X', 0}, {'Z', 0}, {'Y', -2}, {'W', -3}, {'T', 1}} },
{'V', {{'T', 0}, {'D', -3}, {'A', 0}, {'E', -2}, {'I', 3}, {'S', -2}, {'M', 1}, {'Q', -2}, {'L', 1}, {'F', -1}, {'N', -3}, {'R', -3}, {'V', 4}, {'Y', -1}, {'C', -1}, {'G', -3}, {'K', -2}, {'W', -3}, {'H', -3}, {'P', -2}, {'Z', -2}, {'B', -3}, {'X', -1}} },
{'Q', {{'Q', 5}, {'A', -1}, {'R', 1}, {'C', -3}, {'N', 0}, {'D', 0}, {'Y', -1}, {'V', -2}, {'F', -3}, {'L', -2}, {'W', -2}, {'G', -2}, {'M', 0}, {'T', -1}, {'E', 2}, {'K', 1}, {'B', 0}, {'Z', 3}, {'S', 0}, {'H', 0}, {'P', -1}, {'I', -3}, {'X', -1}} },
{'N', {{'A', -2}, {'N', 6}, {'R', 0}, {'L', -3}, {'G', 0}, {'T', 0}, {'M', -2}, {'D', 1}, {'Q', 0}, {'Y', -2}, {'V', -3}, {'F', -3}, {'H', 1}, {'S', 1}, {'C', -3}, {'I', -3}, {'P', -2}, {'X', -1}, {'W', -4}, {'E', 0}, {'K', 0}, {'Z', 0}, {'B', 3}} },
{'Z', {{'Y', -2}, {'Z', 4}, {'P', -1}, {'G', -2}, {'C', -3}, {'K', 1}, {'W', -3}, {'I', -3}, {'V', -2}, {'D', 1}, {'L', -3}, {'H', 0}, {'S', 0}, {'E', 4}, {'A', -1}, {'M', -1}, {'Q', 3}, {'F', -3}, {'T', -1}, {'B', 1}, {'N', 0}, {'R', 0}, {'X', -1}} },
{'H', {{'H', 8}, {'D', -1}, {'G', -2}, {'C', -3}, {'R', 0}, {'N', 1}, {'Q', 0}, {'E', 0}, {'A', -2}, {'S', -1}, {'I', -3}, {'P', -2}, {'X', -1}, {'K', -1}, {'Z', 0}, {'B', 0}, {'W', -2}, {'L', -3}, {'T', -2}, {'M', -2}, {'Y', 2}, {'V', -3}, {'F', -1}} },
{'Y', {{'M', -1}, {'I', -1}, {'E', -2}, {'A', -2}, {'Y', 7}, {'Q', -1}, {'N', -2}, {'F', 3}, {'R', -2}, {'K', -2}, {'G', -3}, {'C', -2}, {'W', 2}, {'S', -2}, {'L', -1}, {'H', 2}, {'D', -3}, {'T', -2}, {'P', -3}, {'Z', -2}, {'B', -3}, {'V', -1}, {'X', -1}} },
{'G', {{'R', -2}, {'N', 0}, {'Q', -2}, {'E', -2}, {'A', 0}, {'D', -1}, {'G', 6}, {'C', -3}, {'K', -2}, {'Z', -2}, {'B', -1}, {'S', 0}, {'H', -2}, {'P', -2}, {'I', -4}, {'X', -1}, {'Y', -3}, {'V', -3}, {'F', -3}, {'W', -2}, {'L', -4}, {'T', -2}, {'M', -3}} },
{'B', {{'Y', -3}, {'S', 0}, {'W', -4}, {'K', 0}, {'G', -1}, {'C', -3}, {'P', -2}, {'L', -4}, {'H', 0}, {'D', 4}, {'I', -3}, {'T', -1}, {'V', -3}, {'Q', 0}, {'M', -3}, {'E', 1}, {'A', -2}, {'R', -1}, {'N', 3}, {'F', -3}, {'B', 4}, {'X', -1}, {'Z', 1}} },
{'E', {{'C', -4}, {'D', 2}, {'Q', 2}, {'A', -1}, {'E', 5}, {'R', 0}, {'N', 0}, {'Y', -2}, {'V', -2}, {'F', -3}, {'W', -3}, {'L', -3}, {'G', -2}, {'T', -1}, {'M', -2}, {'K', 1}, {'Z', 4}, {'B', 1}, {'S', 0}, {'H', 0}, {'P', -1}, {'I', -3}, {'X', -1}} },
{'C', {{'C', 9}, {'A', 0}, {'D', -3}, {'N', -3}, {'R', -3}, {'E', -4}, {'Z', -3}, {'B', -3}, {'Q', -3}, {'S', -1}, {'H', -3}, {'P', -3}, {'I', -1}, {'X', -2}, {'Y', -2}, {'V', -1}, {'K', -3}, {'F', -2}, {'W', -2}, {'L', -1}, {'G', -3}, {'T', -1}, {'M', -1}} },
{'M', {{'R', -1}, {'N', -2}, {'Q', 0}, {'A', -1}, {'E', -2}, {'I', 1}, {'M', 5}, {'D', -3}, {'H', -2}, {'L', 2}, {'C', -1}, {'G', -3}, {'K', -1}, {'Y', -1}, {'V', 1}, {'T', -1}, {'F', 0}, {'W', -1}, {'Z', -1}, {'B', -3}, {'S', -1}, {'P', -2}, {'X', -1}} },
{'T', {{'N', 0}, {'F', -2}, {'R', -1}, {'M', -1}, {'I', -1}, {'A', 0}, {'P', -1}, {'E', -1}, {'Q', -1}, {'T', 5}, {'H', -2}, {'L', -1}, {'D', -1}, {'K', -1}, {'C', -1}, {'G', -2}, {'S', 1}, {'V', 0}, {'X', 0}, {'W', -2}, {'B', -1}, {'Y', -2}, {'Z', -1}} },
{'P', {{'P', 7}, {'D', -1}, {'L', -3}, {'H', -2}, {'G', -2}, {'C', -3}, {'K', -1}, {'R', -2}, {'F', -4}, {'N', -2}, {'Q', -1}, {'E', -1}, {'A', -1}, {'M', -2}, {'I', -3}, {'S', -1}, {'Z', -1}, {'T', -1}, {'B', -2}, {'W', -4}, {'X', -2}, {'Y', -3}, {'V', -2}} },
{'K', {{'K', 5}, {'G', -2}, {'R', 2}, {'H', -1}, {'L', -2}, {'D', -1}, {'Q', 1}, {'I', -3}, {'A', -1}, {'E', 1}, {'C', -3}, {'N', 0}, {'Z', 1}, {'B', 0}, {'S', 0}, {'P', -1}, {'X', -1}, {'Y', -2}, {'V', -2}, {'F', -3}, {'W', -3}, {'T', -1}, {'M', -1}} },
{'I', {{'H', -3}, {'D', -3}, {'G', -4}, {'C', -1}, {'R', -3}, {'N', -3}, {'Q', -3}, {'I', 4}, {'E', -3}, {'A', -1}, {'Y', -1}, {'V', 3}, {'F', 0}, {'Z', -3}, {'W', -3}, {'T', -1}, {'M', 1}, {'B', -3}, {'K', -3}, {'S', -2}, {'P', -3}, {'L', 2}, {'X', -1}} },
{'F', {{'Q', -3}, {'A', -2}, {'E', -3}, {'I', 0}, {'M', 0}, {'R', -3}, {'F', 6}, {'N', -3}, {'C', -2}, {'G', -3}, {'K', -3}, {'D', -3}, {'H', -1}, {'L', 0}, {'W', 1}, {'T', -2}, {'Y', 3}, {'V', -1}, {'S', -2}, {'P', -4}, {'X', -1}, {'Z', -3}, {'B', -3}} },
{'X', {{'L', -1}, {'H', -1}, {'D', -1}, {'X', -1}, {'T', 0}, {'K', -1}, {'G', -1}, {'C', -2}, {'W', -2}, {'S', 0}, {'N', -1}, {'F', -1}, {'B', -1}, {'Z', -1}, {'V', -1}, {'R', -1}, {'P', -2}, {'M', -1}, {'I', -1}, {'E', -1}, {'A', 0}, {'Y', -1}, {'Q', -1}} },
{'D', {{'R', -2}, {'N', 1}, {'A', -2}, {'D', 6}, {'S', 0}, {'H', -1}, {'V', -3}, {'P', -1}, {'I', -3}, {'X', -1}, {'E', 2}, {'C', -3}, {'K', -1}, {'Z', 1}, {'Q', 0}, {'B', 4}, {'W', -4}, {'L', -4}, {'G', -1}, {'T', -1}, {'M', -3}, {'Y', -3}, {'F', -3}} },
{'R', {{'A', -1}, {'R', 5}, {'L', -2}, {'W', -3}, {'G', -2}, {'M', -1}, {'T', -1}, {'D', -2}, {'Q', 1}, {'N', 0}, {'K', 2}, {'Y', -2}, {'V', -3}, {'F', -3}, {'S', -1}, {'H', 0}, {'P', -2}, {'I', -3}, {'X', -1}, {'E', 0}, {'C', -3}, {'B', -1}, {'Z', 0}} },
{'A', {{'A', 4}, {'N', -2}, {'Q', -1}, {'W', -3}, {'Y', -2}, {'V', 0}, {'F', -2}, {'L', -1}, {'C', 0}, {'G', 0}, {'T', 0}, {'M', -1}, {'D', -2}, {'E', -1}, {'R', -1}, {'K', -1}, {'Z', -1}, {'B', -2}, {'S', 1}, {'H', -2}, {'P', -1}, {'I', -1}, {'X', 0}} } };

    //cout << "After blosum init";

    int row, left_cell, current_cell, above_cell, diagonal_cell, match, deletion, insertion, maximum, i, j, max_score;
    vector<int> dp_table( (read_len+1) * (seq_len+1), 0);
    vector<int> traceback_table( (read_len+1) * (seq_len+1), 1);

    vector<int> max_score_coordinates;

    max_score = 0;

    if (penalize_end_gaps) {
        // for needleman-wunsch init the first row and column
        for (int row = 1; row < read_len+1; row++){
            i = (seq_len + 1) * row + 1; // first entry of row
            dp_table[i] = gap_open + (row-1)*gap_score;
        }
        for (int col = 1; col < seq_len+1; col++){
            dp_table[col] = gap_open + (col-1)*gap_score;
        }
    }

    maximum = max_score;

    // cout << "After init" << endl;

    for (int r = 0; r < read_len; r++){
        i = r + 1;
        row = i * (seq_len + 1);
        for (int c = 0; c < seq_len; c++){
            j = c + 1;

            current_cell = row + j;
            left_cell = current_cell - 1;
            above_cell = current_cell - (seq_len + 1);
            diagonal_cell = above_cell - 1;

            // match or missmatch
            if (nuc) {
                if (read[i-1] == seq[j - 1]){
                    match =  dp_table[diagonal_cell] + match_score;
                } else {
                    match = dp_table[diagonal_cell] + mismatch_score;
                }
            } else {
                // the python wrapper multiplies all scores by 2, so here the blosum scores also has to be multiplied
                match =  dp_table[diagonal_cell] + 2*BLOSUM62[read[i-1]][seq[j - 1]];
            }

            if ((i == read_len) && !penalize_end_gaps) {
                deletion = dp_table[left_cell];
            } else if (traceback_table[left_cell] == 0) {
                deletion = dp_table[left_cell] + gap_score;
            } else {
                deletion = dp_table[left_cell] + gap_open;
            }

            if ((j == seq_len) && !penalize_end_gaps) {
                insertion = dp_table[above_cell];
            } else if (traceback_table[above_cell] == 2) {
                insertion = dp_table[above_cell] + gap_score;
            } else {
                insertion = dp_table[above_cell] + gap_open;
            }
            maximum = match;
            traceback_table[current_cell] = 1;

            if (maximum <= deletion){
                maximum = deletion;
                traceback_table[current_cell] = 0;
            }
            if (maximum <= insertion){
                maximum = insertion;
                traceback_table[current_cell] = 2;
            }
            if (sw) {
                if (0 > maximum){
                    maximum = 0;
                }
            }


            // this can be removed and only the the first if kept
            // and only max_score_coordinate = current_cell is required to keep the coordinates for that max score
            if (max_score < maximum){
                max_score = maximum;
                max_score_coordinates.clear();
                max_score_coordinates.push_back(current_cell);
            } else if (max_score == maximum && max_score_coordinates.size() < max_alignments){
                max_score_coordinates.push_back(current_cell);
            }

            dp_table[current_cell] = maximum;
        }
    }

    // Traceback
    // The first loop can be removed and only do the while loop for one max score and one coorindates

    bool continue_traceback;

    if (!sw) {
        max_score_coordinates.clear();
        max_score_coordinates.push_back(current_cell);
    }

    //for (int v: max_score_coordinates){
    //    cout << v << " " << dp_table[v] <<endl;
    //}

    vector<pair <string, string>> alignment_tuples;

    // cout << "Before traceback " << alignment_tuples.size() << endl;

    for (int coord: max_score_coordinates){
        max_score = dp_table[coord];
        int i, j;
        string out_read, out_seq;
        i = coord / (seq_len + 1);
        j = coord % (seq_len + 1);

        if (sw) {
            continue_traceback = (((i != 0) && (j != 0)) && (dp_table[coord] != 0));
        } else {
            continue_traceback = ((i != 0) || (j != 0));
        }

        while (continue_traceback){

            int current_cell = coord;
            int left_cell = current_cell - 1;
            int above_cell = current_cell - (seq_len + 1);
            int diagonal_cell = above_cell - 1;

            if (j == 0) {
                coord = above_cell;
                i = coord / (seq_len + 1);
                j = coord % (seq_len + 1);
                out_read += read[i];
                out_seq += "-";
            } else if (i == 0) {
                coord = left_cell;
                i = coord / (seq_len + 1);
                j = coord % (seq_len + 1);
                out_read += "-";
                out_seq += seq[j];
            } else {

                if (traceback_table[current_cell] == 1){
                    max_score = dp_table[diagonal_cell];
                    coord = diagonal_cell;
                    i = coord / (seq_len + 1);
                    j = coord % (seq_len + 1);
                    out_read += read[i];
                    out_seq += seq[j];

                } else if (traceback_table[current_cell] == 0){
                    max_score = dp_table[left_cell];
                    coord = left_cell;
                    i = coord / (seq_len + 1);
                    j = coord % (seq_len + 1);
                    out_read += "-";
                    out_seq += seq[j];

                } else if (traceback_table[current_cell] == 2){
                    max_score = dp_table[above_cell];
                    coord = above_cell;
                    i = coord / (seq_len + 1);
                    j = coord % (seq_len + 1);
                    out_read += read[i];
                    out_seq += "-";
                } else {
                    cout << "Traceback error";
                    exit(1);
                }
            }
            if (sw) {
                continue_traceback = (((i != 0) && (j != 0)) && (dp_table[coord] != 0));
            } else {
                continue_traceback = ((i != 0) || (j != 0));
            }
        }

        reverse(out_read.begin(), out_read.end());
        reverse(out_seq.begin(), out_seq.end());

        alignment_tuples.push_back({out_read,out_seq});

    }
    return alignment_tuples;
}
