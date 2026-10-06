// digit.cpp — C++ translation of MOPAC 2016 "digit.F90".

#include "digit.h"

double digit(const std::string& s, int istart) {
    const char i0 = '0', i9 = '9', ineg = '-', ipos = '+', idot = '.', ispc = ' ';
    double c1 = 0.0, c2 = 0.0;
    bool sign = true;
    int l = (int)s.size();
    int i = istart - 1;
    int idig = 0;
    for (; i < l; ++i) {
        char n = s[i];
        if (n >= i0 && n <= i9) { idig++; c1 = c1 * 10.0 + (n - i0); }
        else if (n == ineg || n == ipos || n == ispc) { if (n == ineg) sign = false; }
        else if (n == idot) { ++i; break; }
        else { break; }
    }
    double deciml = 1.0;
    for (int j = i; j < l; ++j) {
        char n = s[j];
        if (n >= i0 && n <= i9) { deciml /= 10.0; c2 += (n - i0) * deciml; }
        else if (n != ispc) { break; }
    }
    double v = c1 + c2;
    return sign ? v : -v;
}
