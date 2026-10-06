// geo_diff.cpp — C++ translation of readmo.F90 "geo_diff" (MOPAC 2016).
#include "geo_diff.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

#include "common_arrays_C.h"
#include "elemts_C.h"
#include "molkst_C.h"

using namespace common_arrays_C;
using namespace elemts_C;
using namespace molkst_C;

void geo_diff(double& sum, double& rms, bool prt) {
    static bool first = true;
    std::vector<double> move(numat + 1);
    sum = 0.0;
    rms = 0.0;
    int i, j, k;
    double dum1, sum1;
    for (i = 1; i <= numat; ++i) {
        dum1 = (geo[1][i] - geoa[1][i]) * (geo[1][i] - geoa[1][i]) +
               (geo[2][i] - geoa[2][i]) * (geo[2][i] - geoa[2][i]) +
               (geo[3][i] - geoa[3][i]) * (geo[3][i] - geoa[3][i]);
        rms += dum1;
        dum1 = std::sqrt(dum1);
        move[i] = dum1;
        sum += dum1;
    }
    first = true;
    sum1 = 0.0;
    for (i = 1; i <= std::min(30, numat); ++i) {
        dum1 = -1.0;
        for (j = 1; j <= numat; ++j) {
            if (dum1 < move[j]) {
                dum1 = move[j];
                k = j;
            }
        }
        if (dum1 > 0.1) {
            if (prt && first) {
                first = false;
                std::fprintf(stdout, "\n%22s\n\n%2s\n", "Atoms that move a lot",
                             "Atom No.           Atom Label             GEO_REF Coordinates       Movement    Integral");
            }
            if (prt) {
                sum1 += dum1;
                if (txtatm[k] != " ") {
                    std::fprintf(stdout, "%2d  %4s %8.3f%8.3f%8.3f%12.2f%12.2f\n", k,
                                 ("(" + txtatm[k] + ")").c_str(), geoa[1][k],
                                 geoa[2][k], geoa[3][k], dum1, sum1);
                } else {
                    std::fprintf(stdout, "%15d%4s%23.2f%12.2f\n", k,
                                 elemnt[nat[k]].c_str(), dum1, sum1);
                }
                move[k] = -2.0;
            }
        } else {
            break;
        }
    }
}
