// nuchar.cpp — C++ translation of MOPAC 2016 "nuchar.F90".
// Determines and returns the real values of all numbers found in 'line'.
#include "nuchar.h"
#include "reada.h"

#include <string>

void nuchar(char* line, int l_line, double* value, int& nvalue) {
    int istart[41];
    int i;
    bool leadsp;
    char tab, comma, space;
    tab = char(9);
    comma = ',';
    space = ' ';
    // CLEAN OUT TABS AND COMMAS
    for (i = 1; i <= l_line; ++i) {
        if (line[i - 1] != tab && line[i - 1] != comma) continue;
        line[i - 1] = space;
    }
    // FIND INITIAL DIGIT OF ALL NUMBERS, CHECK FOR LEADING SPACES
    // FOLLOWED BY A CHARACTER
    leadsp = true;
    nvalue = 0;
    for (i = 1; i <= l_line; ++i) {
        if (leadsp && line[i - 1] != space) {
            ++nvalue;
            istart[nvalue] = i;
        }
        leadsp = (line[i - 1] == space);
    }
    // FILL NUMBER ARRAY
    for (i = 1; i <= nvalue; ++i)
        value[i] = reada(std::string(line), istart[i]);
}
