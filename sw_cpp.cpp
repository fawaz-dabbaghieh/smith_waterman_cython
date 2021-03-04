#include "sw_cpp.h"


void smith_waterman(string read, string seq){
    int gap_score = -2;
    int match_score = 1;

    size_t read_len = read.length();
    size_t seq_len = seq.length();


//    int row, left_cell, current_cell, above_cell, diagonal_cell, match, deletion, insertion, maximum, i, j;
    vector<int> dp_table( (read_len+1) * (seq_len+1), 0);
    int current_cell, max_score = 0;

    // you can turn max_score_coorindate to an int to store only one max instead all max scores in dp table
    vector<int> max_score_coordinates;

    int i, j;
    for (int r = 0; r < read_len; r++){
        i = r + 1;
        int row = i * (seq_len + 1);
        for (int c = 0; c < seq_len; c++){
            j = c + 1;

            int maximum = 0;
            int current_cell = row + j;
            int left_cell = current_cell - 1;
            int above_cell = current_cell - seq_len - 1;
            int diagonal_cell = above_cell - 1;
            int match, insertion, deletion;


            // match or missmatch
            if (read[i-1] == seq[j - 1]){
                match =  dp_table[diagonal_cell] + match_score;
            } else {
                match = dp_table[diagonal_cell] - match_score;
            }

            deletion = dp_table[left_cell] + gap_score;
            insertion = dp_table[above_cell] + gap_score;

            if (match > maximum){
                maximum = match;
            }
            if (maximum < deletion){
                maximum = deletion;
            }
            if (maximum < insertion){
                maximum = insertion;
            }
            if (0 > maximum){
                maximum = 0;
            }


            // this can be removed and only the the first if kept
            // and only max_score_coordinate = current_cell is required to keep the coordinates for that max score
            if (max_score < maximum){
                max_score = maximum;
                max_score_coordinates.clear();
                max_score_coordinates.push_back(current_cell);

            } else if (max_score == maximum){
                max_score_coordinates.push_back(current_cell);
            }

            dp_table[current_cell] = maximum;
        }
    }

    // for (int v: max_score_coordinates){
    //     cout << v << " " << dp_table[v] <<endl;
    // }
    // Traceback
    // The first loop can be removed and only do the while loop for one max score and one coorindates
    for (int coord: max_score_coordinates){
        max_score = dp_table[coord];
        int i, j;
        string out_read, out_seq;
        i = coord / (seq_len + 1);
        j = coord % (seq_len + 1);

        while ( ((i != 0) && (j != 0)) && (dp_table[coord] != 0) ){

            int current_cell = coord;
            int left_cell = current_cell - 1;
            int above_cell = current_cell - seq_len - 1;
            int diagonal_cell = above_cell - 1;

            if ((max_score == dp_table[diagonal_cell] + match_score) || (max_score == dp_table[diagonal_cell] - match_score)){
                max_score = dp_table[diagonal_cell];
                coord = diagonal_cell;
                i = coord / (seq_len + 1);
                j = coord % (seq_len + 1);
                out_read += read[i];
                out_seq += seq[j];

            } else if (max_score == dp_table[left_cell] + gap_score){
                max_score = dp_table[left_cell];
                coord = left_cell;
                i = coord / (seq_len + 1);
                j = coord % (seq_len + 1);
                out_read += '-';
                out_seq += seq[j];

            } else if (max_score == dp_table[above_cell] + gap_score){
                max_score = dp_table[above_cell];
                coord = above_cell;
                i = coord / (seq_len + 1);
                j = coord % (seq_len + 1);
                out_read += read[i];
                out_seq += '-';
            }
        }

        cout << endl;
        for (int i=out_read.length() -1; i>=0; i--){
            cout << out_read[i];
        }
        cout << endl;
        for (int i=out_seq.length() -1; i>=0; i--){
            cout << out_seq[i];
        }
        // reverse(out_read.begin(), out_read.end()); 
        // reverse(out_seq.begin(), out_seq.end());
        // cout << out_read << endl;
        // cout << out_seq <<  endl;
    }
}
