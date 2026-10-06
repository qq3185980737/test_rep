// mult33.cpp
#include "mult33.h"
#include "symmetry_C.h"

using namespace symmetry_C;

void mult33(const double* fmat, int iplace) {
    // fmat is Fortran 3x3 column-major: fmat(i,l) = fmat[(l-1)*3 + (i-1)].
    double help[4][4] = {};
    for (int i = 1; i <= 3; ++i)
        for (int j = 1; j <= 3; ++j) {
            help[i][j] = 0.0;
            for (int k = 1; k <= 3; ++k) {
                double s = 0.0;
                for (int l = 1; l <= 3; ++l)
                    s += fmat[(l-1)*3 + (i-1)] * fmat[(k-1)*3 + (j-1)] * elem[l][k][iplace];
                help[i][j] += s;
            }
        }
    for (int i = 1; i <= 3; ++i)
        for (int j = 1; j <= 3; ++j)
            elem[i][j][iplace] = help[i][j];
}
