// dtrans.cpp — C++ translation of MOPAC 2016 "dtrans.F90".

#include "dtrans.h"

#include "dtran2.h"
#include "symmetry_C.h"

using namespace symmetry_C;

void dtrans(std::vector<double>& d, int ioper, bool& first,
            const std::vector<std::vector<double>>& orient) {
    static std::vector<std::vector<std::vector<double>>> t1(13,
        std::vector<std::vector<double>>(5, std::vector<double>(5, 0.0)));
    if (first) {
        first = false;
        std::vector<std::vector<double>> s(3, std::vector<double>(3, 0.0));
        for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) s[i][j] = orient[i][j];
        dtran2(s, t1, 1);
        for (int k = 2; k <= nclass; ++k) {
            for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) s[i][j] = elem[i+1][j+1][k];
            dtran2(s, t1, k);
        }
    }
    std::vector<double> h(6, 0.0);
    for (int i = 1; i <= 5; ++i) {
        double s = 0.0;
        for (int j = 1; j <= 5; ++j) s += t1[ioper][i-1][j-1] * d[j];
        h[i] = s;
    }
    for (int i = 1; i <= 5; ++i) {
        double s = 0.0;
        for (int j = 1; j <= 5; ++j) s += t1[ioper][j-1][i-1] * h[j];
        d[i] = s;
    }
}
