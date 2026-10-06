// prtgra.cpp
#include "prtgra.h"
#include <vector>
#include <cmath>
#include <cstdio>
#include "common_arrays_C.h"
#include "molkst_C.h"
using common_arrays_C::grad;
using common_arrays_C::xparam;
using common_arrays_C::txtatm;
using common_arrays_C::labels;
using common_arrays_C::loc;
using common_arrays_C::na;
using molkst_C::natoms;
using molkst_C::nvar;
using molkst_C::maxtxt;
using molkst_C::na1;
using molkst_C::keywrd;
extern const char* elemnt(int);
extern int iw;

void prtgra() {
    double glim = 0.1;
    auto idx = [&](const std::string& s, const char* q){ return s.find(q); };
    int i = (int)(idx(keywrd," GRAD=") != std::string::npos) +
            (int)(idx(keywrd,"DERIV") != std::string::npos) +
            (int)(idx(keywrd," GRAD") != std::string::npos);
    if (i != 0) glim = -0.1;
    std::printf("\n");
    int l = 1;
    double tmp[4][3];
    for (int iatom = 1; iatom <= natoms; ++iatom) {
        int jcnt = 0;
        double gmax = 0;
        if (l <= nvar) {
            for (int m = 1; m <= 3; ++m) {
                if (loc[1][l] == iatom && loc[2][l] == m) {
                    ++jcnt;
                    tmp[jcnt][1] = grad[l];
                    tmp[jcnt][2] = xparam[l];
                    if (na1 != 99 && m > 1 && na[iatom] != 0)
                        tmp[jcnt][2] *= 180.0 / 3.14159265358979;
                    if (std::fabs(grad[l]) > gmax) gmax = std::fabs(grad[l]);
                    ++l;
                }
            }
        }
        if (jcnt != 0 && glim < gmax) {
            double sum = 0;
            if (jcnt == 3 && na[iatom] == 0)
                sum = std::sqrt(tmp[1][1]*tmp[1][1]+tmp[2][1]*tmp[2][1]+tmp[3][1]*tmp[3][1]);
            std::printf("%5d %s", iatom, elemnt(labels[iatom]));
            for (int n=1;n<=2;++n) for (int m=1;m<=jcnt;++m)
                std::printf(" %9.3f", tmp[m][n]);
            if (jcnt == 3 && na[iatom] == 0) std::printf(" %12.3f", sum);
            std::printf("\n");
        }
    }
    std::printf("\n");
}
