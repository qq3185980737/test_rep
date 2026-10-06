// reada.cpp — C++ translation of "reada.F90".
// Extract a real number from a string starting at a 1-based position.

#include "reada.h"
#include "digit.h"

#include <cmath>
#include <string>

double reada(const std::string& string, int istart) {
    const int i0 = '0', i9 = '9', idot = '.', ineg = '-', ipos = '+';
    const int icapd = 'D', icape = 'E', ismld = 'd', ismle = 'e';
    int l = (int)string.size();

    // Find start of numeric field (Fortran istart is 1-based).
    int i = -1;
    for (int ii = istart; ii <= l; ++ii) {
        int n = (unsigned char)string[ii - 1];
        if (n >= i0 && n <= i9) { i = ii; break; }
        if (n == ineg || n == ipos) {
            if (ii + 1 > l) return 0.0;
            int nn = (unsigned char)string[ii];  // C++ index for Fortran ii+1
            if (nn >= i0 && nn <= i9) { i = ii; break; }
        }
        if (n != idot) continue;
        if (ii + 1 > l) return 0.0;
        int nn = (unsigned char)string[ii];
        if (nn >= i0 && nn <= i9) { i = ii; break; }
    }
    if (i < 0) return 0.0;

    // Find end of numeric field.
    bool expnnt = false;
    int j = l + 1;
    for (int jj = i + 1; jj <= l; ++jj) {
        int n = (unsigned char)string[jj - 1];
        if (n >= i0 && n <= i9) continue;
        if (n == ineg || n == ipos) {
            if (!expnnt) { j = jj; break; }
            if (jj + 1 > l) { j = jj; break; }
            int nn = (unsigned char)string[jj];
            if (nn >= i0 && nn <= i9) continue;
        }
        if (n == idot) {
            if (jj + 1 > l) { j = jj; break; }
            int nn = (unsigned char)string[jj];
            if (nn >= i0 && nn <= i9) continue;
            if (nn == icape || nn == ismle || nn == icapd || nn == ismld) continue;
        }
        if (n == icape || n == ismle || n == icapd || n == ismld) {
            if (expnnt) { j = jj; break; }
            expnnt = true;
            continue;
        }
        j = jj;
        break;
    }
    // Fortran: char at j-1.
    if (j - 2 >= 0 && j - 2 < l) {
        int last = (unsigned char)string[j - 2];
        if (last == icape || last == ismle || last == icapd || last == ismld) --j;
    }

    // Numeric field is Fortran string[i..j-1].
    std::string field = string.substr(i - 1, (j - 1) - i + 1);
    size_t ep = std::string::npos;
    for (size_t k = 0; k < field.size(); ++k) {
        char c = field[k];
        if (c == 'e' || c == 'E' || c == 'd' || c == 'D') { ep = k; break; }
    }
    if (ep == std::string::npos) return digit(field, 1);
    double mantissa = digit(field.substr(0, ep), 1);
    double exponent = digit(field.substr(ep + 1), 1);
    return mantissa * std::pow(10.0, exponent);
}
