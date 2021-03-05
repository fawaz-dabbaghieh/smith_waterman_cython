#ifndef SW_CPP_H
#define SW_CPP_H

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <algorithm>

using namespace std;


vector<pair<string, string>> smith_waterman(std::string read, std::string seq, int gap_open, int gap_score, bool penalize_end_gaps,
                    bool sw, bool nuc, int nuc_match, int nuc_mismatch, int max_alignments);


# endif // SW_CPP_H
