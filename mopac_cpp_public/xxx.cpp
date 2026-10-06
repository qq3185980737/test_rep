// xxx.cpp
#include "xxx.h"
void xxx(char type, int i, int j, int k, int l, std::string& r) {
    // Fortran r*13: r(1)=type, digits follow. C++ 0-based: r[0]=type.
    r.clear();
    r.push_back(type);
    int ijk[4] = {i, j, k, l};
    for (int loop = 0; loop < 4; ++loop) {
        int ii = ijk[loop];
        if (ii == 0) continue;
        int i2 = ii / 10;
        if (i2 != 0) { r.push_back(char('0' + i2)); ii -= i2 * 10; }
        r.push_back(char('0' + ii));
    }
}
