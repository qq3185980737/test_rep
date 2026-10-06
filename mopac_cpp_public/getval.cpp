// getval.cpp — C++ translation of MOPAC 2016 "getval.F90".
// F90: numeric iff ch1 AND ch2 are both NOT uppercase 'A'..'Z'; x = reada(line,1).
#include "getval.h"
#include "reada.h"
#include <string>

void getval(const std::string& line, double& x, std::string& t) {
    char ch1 = line.size() > 0 ? line[0] : ' ';
    char ch2 = line.size() > 1 ? line[1] : ' ';
    bool letter1 = (ch1 >= 'A' && ch1 <= 'Z');
    bool letter2 = (ch2 >= 'A' && ch2 <= 'Z');
    if (!letter1 && !letter2) {
        x = reada(line, 1);
        t = " ";
    } else {
        size_t sp = line.find(' ');
        t = line.substr(0, sp);
        x = -999.0;
    }
}
